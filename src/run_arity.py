"""
run_arity.py — compute per-protein arity signatures from the v4 tarball.

Streams the tarball, extracts Cα coordinates via gemmi, computes
(arity_25, arity_75) for each protein, writes a parquet with:
    uniprot_id, filename, fragment, n_residues, arity_25, arity_75

This file is then merged with proteome_full.parquet on (uniprot_id, filename).
"""
from __future__ import annotations

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

_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from arity import arity_signature  # noqa: E402

log = logging.getLogger(__name__)

_ROOT     = Path(__file__).resolve().parent.parent
_TAR_PATH = _ROOT / "data" / "hsapiens" / "UP000005640_9606_HUMAN_v4.tar"
_OUT_PATH = _ROOT / "data" / "results" / "proteome_full_arity.parquet"
_FAIL_LOG = _ROOT / "data" / "results" / "proteome_arity_failures.log"

_AF_RE = re.compile(r"AF-(?P<uid>[^-]+)-F(?P<frag>\d+)-model_v\d+", re.I)

_SCHEMA = pa.schema([
    ("uniprot_id", pa.string()),
    ("filename",   pa.string()),
    ("fragment",   pa.int16()),
    ("n_residues", pa.int32()),
    ("arity_25",   pa.int16()),
    ("arity_75",   pa.int16()),
])


def _worker(args: tuple) -> tuple[dict | None, str]:
    content_bytes, filename = args
    try:
        pdb_str = content_bytes.decode("utf-8", errors="replace")
    except Exception as exc:
        return None, f"{filename}: decode error: {exc}"

    try:
        import gemmi
        st    = gemmi.read_pdb_string(pdb_str)
        model = st[0]
    except Exception as exc:
        return None, f"{filename}: gemmi error: {exc}"

    coords: list[list[float]] = []
    for chain in model:
        for res in chain:
            ca = res.find_atom("CA", "\0")
            if ca is not None:
                pos = ca.pos
                coords.append([pos.x, pos.y, pos.z])

    n = len(coords)
    if n < 2:
        return None, f"{filename}: only {n} Cα atoms"

    ca_coords = np.array(coords, dtype=np.float64)
    a25, a75  = arity_signature(ca_coords)

    m    = _AF_RE.search(filename)
    uid  = m.group("uid")  if m else Path(filename).stem
    frag = int(m.group("frag")) if m else 1

    return {
        "uniprot_id": uid,
        "filename":   Path(filename).name,
        "fragment":   frag,
        "n_residues": n,
        "arity_25":   a25,
        "arity_75":   a75,
    }, ""


def run(tar_path: Path = _TAR_PATH,
        out_path: Path = _OUT_PATH,
        workers: int = 6,
        chunk_size: int = 200,
        flush_every: int = 1000) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    n_failed  = 0
    pending: list[dict] = []

    with (
        pq.ParquetWriter(str(out_path), _SCHEMA) as writer,
        open(_FAIL_LOG, "w") as fail_fh,
        tarfile.open(str(tar_path), mode="r|") as tar,
        multiprocessing.Pool(processes=workers) as pool,
    ):
        pbar  = tqdm(desc="Arity signatures", unit="pdb", smoothing=0.05)
        chunk: list[tuple] = []

        def _flush_chunk():
            nonlocal n_failed
            if not chunk:
                return
            for row, err in pool.map(_worker, chunk):
                if row is not None:
                    pending.append(row)
                else:
                    n_failed += 1
                    if err:
                        fail_fh.write(err + "\n")
            pbar.update(len(chunk))
            chunk.clear()

        def _flush_parquet():
            nonlocal n_written
            if not pending:
                return
            writer.write_table(pa.table(
                {col: [r[col] for r in pending] for col in _SCHEMA.names},
                schema=_SCHEMA,
            ))
            n_written += len(pending)
            pending.clear()

        for member in tar:
            name = member.name
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
            chunk.append((raw, name))
            if len(chunk) >= chunk_size:
                _flush_chunk()
                if len(pending) >= flush_every:
                    _flush_parquet()

        _flush_chunk()
        _flush_parquet()
        pbar.close()

    log.info("Done. %d written, %d failed.", n_written, n_failed)
    return n_written


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    n = run()
    print(f"\nTotal written: {n}")
