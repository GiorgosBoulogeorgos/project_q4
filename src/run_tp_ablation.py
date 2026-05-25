"""
run_tp_ablation.py — reprocess the proteome tarball computing PLM at three
thresholds (t_p = 0.020, 0.025, 0.030) in a single pass per protein.

Output: data/results/proteome_full_tp_ablation.parquet
Columns: uniprot_id, filename, fragment, n_residues,
         f_cp_plus, H_p, PLM_020, PLM_025, PLM_030
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

from filtration import build_pd_and_ncc   # noqa: E402
from plm import plm_multi                 # noqa: E402
from statistics import f_cp_plus, persistence_entropy  # noqa: E402

log = logging.getLogger(__name__)

_ROOT     = Path(__file__).resolve().parent.parent
_TAR_PATH = _ROOT / "data" / "hsapiens" / "UP000005640_9606_HUMAN_v4.tar"
_OUT_PATH = _ROOT / "data" / "results" / "proteome_full_tp_ablation.parquet"
_FAIL_LOG = _ROOT / "data" / "results" / "proteome_tp_ablation_failures.log"

_AF_RE = re.compile(r"AF-(?P<uid>[^-]+)-F(?P<frag>\d+)-model_v\d+", re.I)

_T_PS = (0.020, 0.025, 0.030)

_SCHEMA = pa.schema([
    ("uniprot_id", pa.string()),
    ("filename",   pa.string()),
    ("fragment",   pa.int16()),
    ("n_residues", pa.int32()),
    ("f_cp_plus",  pa.float32()),
    ("H_p",        pa.float32()),
    ("PLM_020",    pa.int16()),
    ("PLM_025",    pa.int16()),
    ("PLM_030",    pa.int16()),
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

    plddt_vals: list[float] = []
    for chain in model:
        for res in chain:
            ca = res.find_atom("CA", "\0")
            if ca is not None:
                plddt_vals.append(float(ca.b_iso))

    n = len(plddt_vals)
    if n < 2:
        return None, f"{filename}: only {n} Cα atoms"

    plddt = np.array(plddt_vals, dtype=np.float64)

    try:
        pd_arr, ncc_curve = build_pd_and_ncc(plddt)
    except Exception as exc:
        return None, f"{filename}: filtration error: {exc}"

    m    = _AF_RE.search(filename)
    uid  = m.group("uid")  if m else Path(filename).stem
    frag = int(m.group("frag")) if m else 1

    plm_vals = plm_multi(ncc_curve, n_residues=n, t_ps=_T_PS)

    return {
        "uniprot_id": uid,
        "filename":   Path(filename).name,
        "fragment":   frag,
        "n_residues": n,
        "f_cp_plus":  float(f_cp_plus(pd_arr)),
        "H_p":        float(persistence_entropy(pd_arr)),
        "PLM_020":    plm_vals[0],
        "PLM_025":    plm_vals[1],
        "PLM_030":    plm_vals[2],
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
        pbar  = tqdm(desc="t_p ablation", unit="pdb", smoothing=0.05)
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
            table = pa.table(
                {col: [r[col] for r in pending] for col in _SCHEMA.names},
                schema=_SCHEMA,
            )
            writer.write_table(table)
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
