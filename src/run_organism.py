"""
run_organism.py — run the full Q4 pipeline for any AlphaFold-DB model organism.

This is a thin orchestrator over the existing, organism-agnostic pipeline.  The
scientific core (filtration, statistics, PLM, arity) operates on bare pLDDT/Cα
arrays and knows nothing about species, so adding an organism only means:

  1. resolving its UniProt proteome ID and downloading the tarball, and
  2. pointing the three existing entry points at organism-named outputs.

For each organism it produces three parquets that mirror the H. sapiens set
consumed by notebooks/04_proteome_figures7_8*.ipynb:

    data/results/proteome_<org>_<ver>_full.parquet         (main stats + PLM)
    data/results/proteome_<org>_<ver>_tp_ablation.parquet  (PLM at 3 thresholds)
    data/results/proteome_<org>_<ver>_arity.parquet        (Q1 arity signatures)

Usage
-----
    python src/run_organism.py MOUSE [--version v6] [--workers 6] [--t-p 0.025]
    python src/run_organism.py YEAST --tar data/yeast/UP000002311_559292_YEAST_v6.tar

Supported species codes: MOUSE, RAT, YEAST, PSEAE, DROME (and HUMAN, for
completeness). PSEAE is P. aeruginosa PAO1, a bacterial out-of-clade probe;
DROME is D. melanogaster, a non-mammalian (invertebrate) animal probe.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import sys
from pathlib import Path

# ── ensure src/ on path so the imported drivers find filtration/plm/etc. ─────
_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import run_arity                              # noqa: E402  (module: override _FAIL_LOG)
import run_tp_ablation                        # noqa: E402  (module: override _FAIL_LOG)
from download import _PROTEOME_IDS, _resume_get  # noqa: E402
from run_full_proteome import process_tar     # noqa: E402

log = logging.getLogger(__name__)

_ROOT      = Path(__file__).resolve().parent.parent
_RESULTS   = _ROOT / "data" / "results"
_FTP_BASE  = "https://ftp.ebi.ac.uk/pub/databases/alphafold"

# data/<dir> per species, mirroring the existing data/hsapiens convention.
_DATA_DIR = {
    "HUMAN": "hsapiens",
    "MOUSE": "mouse",
    "RAT":   "rat",
    "YEAST": "yeast",
    "PSEAE": "pseae",   # P. aeruginosa PAO1 — bacterial generalisation probe
    "DROME": "drome",   # D. melanogaster — invertebrate (non-mammalian animal) probe
}


def ensure_tarball(species: str, version: str, tar_dir: Path) -> Path:
    """Download the proteome tarball for *species*/*version* if not present."""
    proteome_id = _PROTEOME_IDS[species]
    tar_name    = f"{proteome_id}_{version}.tar"
    dest        = tar_dir / tar_name
    tar_dir.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > 1_000_000:
        log.info("Tarball already present: %s (%.2f GB)",
                 dest, dest.stat().st_size / 1e9)
        return dest

    url = f"{_FTP_BASE}/{version}/{tar_name}"
    log.info("Downloading %s → %s", url, dest)
    _resume_get(url, dest)
    return dest


def run_species(
    species: str,
    version: str = "v6",
    workers: int = 6,
    t_p: float = 0.025,
    chunk_size: int = 200,
    flush_every: int = 1000,
    tar_path: Path | None = None,
) -> dict[str, int]:
    """
    Run the full pipeline (main + t_p ablation + arity) for one organism.

    Returns a dict of {stage: n_proteins_written}.
    """
    species = species.upper()
    if species not in _PROTEOME_IDS:
        raise ValueError(
            f"Unknown species {species!r}. Supported: {sorted(_PROTEOME_IDS)}"
        )

    # output filename token; HUMAN is supported for completeness, but its
    # canonical v6 run is run_v6_proteome.py (-> proteome_v6_full.parquet).
    org      = "hsapiens" if species == "HUMAN" else species.lower()
    tar_dir  = _ROOT / "data" / _DATA_DIR[species]
    tar_path = tar_path or ensure_tarball(species, version, tar_dir)

    _RESULTS.mkdir(parents=True, exist_ok=True)
    stem = f"proteome_{org}_{version}"
    out_full = _RESULTS / f"{stem}_full.parquet"
    out_abl  = _RESULTS / f"{stem}_tp_ablation.parquet"
    out_ari  = _RESULTS / f"{stem}_arity.parquet"

    results: dict[str, int] = {}

    # ── 1. main stats + PLM(t_p) ─────────────────────────────────────────────
    log.info("[%s] stage 1/3 — main stats → %s", org, out_full.name)
    results["full"] = process_tar(
        tar_path    = tar_path,
        out_path    = out_full,
        fail_log    = _RESULTS / f"{stem}_failures.log",
        t_p         = t_p,
        workers     = workers,
        chunk_size  = chunk_size,
        flush_every = flush_every,
    )

    # ── 2. t_p ablation (PLM at 0.020 / 0.025 / 0.030) ───────────────────────
    # run_tp_ablation.run() writes to a module-level _FAIL_LOG; redirect it so
    # organisms don't clobber each other's logs.
    log.info("[%s] stage 2/3 — t_p ablation → %s", org, out_abl.name)
    run_tp_ablation._FAIL_LOG = _RESULTS / f"{stem}_tp_ablation_failures.log"
    results["tp_ablation"] = run_tp_ablation.run(
        tar_path    = tar_path,
        out_path    = out_abl,
        workers     = workers,
        chunk_size  = chunk_size,
        flush_every = flush_every,
    )

    # ── 3. Q1 arity signatures ───────────────────────────────────────────────
    log.info("[%s] stage 3/3 — arity → %s", org, out_ari.name)
    run_arity._FAIL_LOG = _RESULTS / f"{stem}_arity_failures.log"
    results["arity"] = run_arity.run(
        tar_path    = tar_path,
        out_path    = out_ari,
        workers     = workers,
        chunk_size  = chunk_size,
        flush_every = flush_every,
    )

    log.info("[%s] done: %s", org,
             ", ".join(f"{k}={v}" for k, v in results.items()))
    return results


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Run the full Q4 pipeline for an AlphaFold-DB model organism."
    )
    parser.add_argument("species",
                        help="MOUSE | RAT | YEAST | PSEAE | DROME | HUMAN")
    parser.add_argument("--version", default="v6",
                        help="AlphaFold-DB version (default v6)")
    parser.add_argument("--workers", type=int,   default=6)
    parser.add_argument("--t-p",     type=float, default=0.025)
    parser.add_argument("--chunk",   type=int,   default=200)
    parser.add_argument("--flush",   type=int,   default=1000)
    parser.add_argument("--tar",     type=Path,  default=None,
                        help="Path to an existing tarball (skips download)")
    args = parser.parse_args()

    res = run_species(
        species     = args.species,
        version     = args.version,
        workers     = args.workers,
        t_p         = args.t_p,
        chunk_size  = args.chunk,
        flush_every = args.flush,
        tar_path    = args.tar,
    )
    print("\nWritten:", ", ".join(f"{k}={v}" for k, v in res.items()))
