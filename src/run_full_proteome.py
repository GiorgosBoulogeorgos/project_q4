"""
run_full_proteome.py — stream-process the full H. Sapiens AlphaFold-DB v4 tarball.

Usage
-----
    python src/run_full_proteome.py [--workers N] [--t-p T_P] [--chunk CHUNK]

Workflow
--------
1. Download UP000005640_9606_HUMAN_v4.tar (~5 GB, resumed if interrupted).
2. Stream the tar member-by-member; each member is a gzip-compressed PDB.
3. Decompress each member in memory; send bytes to a process-pool worker.
4. Worker: parse with gemmi.read_pdb_string, run build_pd_and_ncc (which by
   default rounds pLDDT to integers in [0, 100] before the filtration —
   the default behaviour as of the integer-pLDDT switch), compute
   n_residues, f⁺_cp, mean_persistence, H_p, ncc_max, PLM(t_p).
5. Flush to data/results/proteome_full.parquet every FLUSH_EVERY proteins
   using PyArrow ParquetWriter (row-group append).
6. Log failures to data/results/proteome_failures.log.

Each AlphaFold-DB fragment file (F1, F2, …) is a separate row, matching
the paper's protein-count semantics.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import multiprocessing
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

# ── ensure src/ on path so workers can import project modules ───────────────
_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from download import _resume_get         # noqa: E402
from filtration import build_pd_and_ncc # noqa: E402
from plm import plm as compute_plm      # noqa: E402
from statistics import (                 # noqa: E402
    f_cp_plus, mean_persistence, ncc_max, persistence_entropy,
)

log = logging.getLogger(__name__)

_ROOT     = Path(__file__).resolve().parent.parent
_TAR_DIR  = _ROOT / "data" / "hsapiens"
_TAR_NAME = "UP000005640_9606_HUMAN_v4.tar"
_TAR_URL  = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/{_TAR_NAME}"
_OUT_PATH = _ROOT / "data" / "results" / "proteome_full.parquet"
_FAIL_LOG = _ROOT / "data" / "results" / "proteome_failures.log"

_AF_RE   = re.compile(r"AF-(?P<uid>[^-]+)-F(?P<frag>\d+)-model_v\d+", re.I)
_THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

_PARQUET_SCHEMA = pa.schema([
    ("uniprot_id",       pa.string()),
    ("filename",         pa.string()),
    ("fragment",         pa.int16()),
    ("n_residues",       pa.int32()),
    ("f_cp_plus",        pa.float32()),
    ("mean_persistence", pa.float32()),
    ("H_p",             pa.float32()),
    ("ncc_max",          pa.int32()),
    ("PLM",              pa.int16()),
])

# ---------------------------------------------------------------------------
# Module-level worker (must be importable by child processes)
# ---------------------------------------------------------------------------

def _compute_from_bytes(args: tuple) -> tuple[dict | None, str]:
    """
    Parse a gzip-compressed AlphaFold PDB from bytes and compute statistics.

    Parameters
    ----------
    args : (content_bytes, filename, t_p)

    Returns
    -------
    (row_dict or None, error_message_or_empty_string)
    """
    content_bytes, filename, t_p = args
    try:
        pdb_str = content_bytes.decode("utf-8", errors="replace")
    except Exception as exc:
        return None, f"{filename}: decode error: {exc}"

    try:
        import gemmi
        st = gemmi.read_pdb_string(pdb_str)
        model = st[0]
    except Exception as exc:
        return None, f"{filename}: gemmi parse error: {exc}"

    plddt_vals: list[float] = []
    for chain in model:
        for res in chain:
            ca = res.find_atom("CA", "\0")
            if ca is None:
                continue
            plddt_vals.append(float(ca.b_iso))

    n = len(plddt_vals)
    if n < 2:
        return None, f"{filename}: only {n} Cα atoms, skipped"

    plddt = np.array(plddt_vals, dtype=np.float64)

    try:
        pd_arr, ncc_curve = build_pd_and_ncc(plddt)
    except Exception as exc:
        return None, f"{filename}: filtration error: {exc}"

    m = _AF_RE.search(filename)
    uid  = m.group("uid")  if m else Path(filename).stem
    frag = int(m.group("frag")) if m else 1

    try:
        plm_val = compute_plm(ncc_curve, n_residues=n, t_p=t_p)
    except Exception:
        plm_val = -1

    row = {
        "uniprot_id":       uid,
        "filename":         Path(filename).name,
        "fragment":         frag,
        "n_residues":       n,
        "f_cp_plus":        float(f_cp_plus(pd_arr)),
        "mean_persistence": float(mean_persistence(pd_arr)),
        "H_p":             float(persistence_entropy(pd_arr)),
        "ncc_max":          int(ncc_max(ncc_curve)),
        "PLM":              int(plm_val),
    }
    return row, ""


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def ensure_tarball(tar_dir: Path = _TAR_DIR) -> Path:
    """Download the proteome tarball if not already present; return its path."""
    tar_dir.mkdir(parents=True, exist_ok=True)
    dest = tar_dir / _TAR_NAME
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log.info("Tarball already present: %s (%.1f GB)",
                 dest, dest.stat().st_size / 1e9)
        return dest
    log.info("Downloading %s → %s", _TAR_URL, dest)
    _resume_get(_TAR_URL, dest)
    return dest


# ---------------------------------------------------------------------------
# Core streaming pipeline
# ---------------------------------------------------------------------------

def process_tar(
    tar_path: Path,
    out_path: Path,
    fail_log: Path = _FAIL_LOG,
    t_p: float = 0.025,
    workers: int = 6,
    chunk_size: int = 200,
    flush_every: int = 1000,
) -> int:
    """
    Stream the proteome tarball and write statistics to a Parquet file.

    Returns the total number of proteins written.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_failed  = 0
    pending_rows: list[dict] = []

    with (
        pq.ParquetWriter(str(out_path), _PARQUET_SCHEMA) as pq_writer,
        open(fail_log, "w") as fail_fh,
        tarfile.open(str(tar_path), mode="r|") as tar,
        multiprocessing.Pool(processes=workers) as pool,
    ):
        pbar = tqdm(desc="Processing proteins", unit="pdb", smoothing=0.05)

        chunk: list[tuple] = []

        def _flush_chunk() -> None:
            nonlocal n_written, n_failed
            if not chunk:
                return
            results = pool.map(_compute_from_bytes, chunk)
            for row, err in results:
                if row is not None:
                    pending_rows.append(row)
                else:
                    n_failed += 1
                    if err:
                        fail_fh.write(err + "\n")
            pbar.update(len(chunk))
            chunk.clear()

        def _flush_parquet() -> None:
            nonlocal n_written
            if not pending_rows:
                return
            table = pa.table(
                {col: [r[col] for r in pending_rows]
                 for col in _PARQUET_SCHEMA.names},
                schema=_PARQUET_SCHEMA,
            )
            pq_writer.write_table(table)
            n_written += len(pending_rows)
            pending_rows.clear()
            log.info("Flushed %d proteins (total %d)", len(table), n_written)

        for member in tar:
            name = member.name
            # AlphaFold-DB v4 tar members: AF-{uid}-F{n}-model_v4.pdb.gz
            if not ("-model_v" in name and
                    (name.endswith(".pdb.gz") or name.endswith(".pdb"))):
                continue

            f = tar.extractfile(member)
            if f is None:
                continue

            raw = f.read()
            if name.endswith(".gz"):
                try:
                    raw = gzip.decompress(raw)
                except Exception as exc:
                    fail_fh.write(f"{name}: gzip error: {exc}\n")
                    n_failed += 1
                    continue

            chunk.append((raw, name, t_p))

            if len(chunk) >= chunk_size:
                _flush_chunk()
                if len(pending_rows) >= flush_every:
                    _flush_parquet()

        # drain remaining
        _flush_chunk()
        _flush_parquet()
        pbar.close()

    log.info("Done. %d proteins written, %d failures.", n_written, n_failed)
    return n_written


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Full H. Sapiens proteome pipeline (streaming tar)"
    )
    parser.add_argument("--workers",   type=int,   default=6,
                        help="Parallel worker processes (default 6)")
    parser.add_argument("--t-p",       type=float, default=0.025,
                        help="PLM threshold (default 0.025)")
    parser.add_argument("--chunk",     type=int,   default=200,
                        help="Proteins per pool.map call (default 200)")
    parser.add_argument("--flush",     type=int,   default=1000,
                        help="Proteins per parquet row-group flush (default 1000)")
    parser.add_argument("--out",       type=Path,  default=_OUT_PATH,
                        help=f"Output parquet (default {_OUT_PATH})")
    parser.add_argument("--tar",       type=Path,  default=None,
                        help="Path to existing tarball (skips download)")
    args = parser.parse_args()

    tar_path = args.tar or ensure_tarball()
    n = process_tar(
        tar_path=tar_path,
        out_path=args.out,
        t_p=args.t_p,
        workers=args.workers,
        chunk_size=args.chunk,
        flush_every=args.flush,
    )
    print(f"\nTotal proteins written: {n}")
