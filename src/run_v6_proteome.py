"""
run_v6_proteome.py — stream-process the full H. Sapiens AlphaFold-DB v6 tarball.

Usage
-----
    python src/run_v6_proteome.py [--workers N] [--t-p T_P] [--chunk CHUNK]

Identical pipeline to run_full_proteome.py; only the tarball URL and output
paths differ (v4 → v6).  The _AF_RE regex already matches any model_vN, so
no parsing changes are needed.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import sys
from pathlib import Path

# ── re-use every symbol from the v4 module except the three constants ────────
_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from run_full_proteome import (   # noqa: E402
    _compute_from_bytes,          # worker (must stay importable by child procs)
    _PARQUET_SCHEMA,
    process_tar,
)
from download import _resume_get  # noqa: E402

log = logging.getLogger(__name__)

_ROOT     = Path(__file__).resolve().parent.parent
_TAR_DIR  = _ROOT / "data" / "hsapiens"
_TAR_NAME = "UP000005640_9606_HUMAN_v6.tar"
_TAR_URL  = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/{_TAR_NAME}"
_OUT_PATH = _ROOT / "data" / "results" / "proteome_v6_full.parquet"
_FAIL_LOG = _ROOT / "data" / "results" / "proteome_v6_failures.log"


def ensure_tarball(tar_dir: Path = _TAR_DIR) -> Path:
    tar_dir.mkdir(parents=True, exist_ok=True)
    dest = tar_dir / _TAR_NAME
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log.info("Tarball already present: %s (%.1f GB)",
                 dest, dest.stat().st_size / 1e9)
        return dest
    log.info("Downloading %s → %s", _TAR_URL, dest)
    _resume_get(_TAR_URL, dest)
    return dest


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Full H. Sapiens proteome pipeline — AlphaFold-DB v6"
    )
    parser.add_argument("--workers", type=int,   default=6)
    parser.add_argument("--t-p",     type=float, default=0.025)
    parser.add_argument("--chunk",   type=int,   default=200)
    parser.add_argument("--flush",   type=int,   default=1000)
    parser.add_argument("--out",     type=Path,  default=_OUT_PATH)
    parser.add_argument("--tar",     type=Path,  default=None,
                        help="Path to existing tarball (skips download)")
    args = parser.parse_args()

    tar_path = args.tar or ensure_tarball()
    n = process_tar(
        tar_path  = tar_path,
        out_path  = args.out,
        fail_log  = _FAIL_LOG,
        t_p       = args.t_p,
        workers   = args.workers,
        chunk_size= args.chunk,
        flush_every=args.flush,
    )
    print(f"\nTotal proteins written: {n}")
