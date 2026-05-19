"""
process_proteome.py — batch pipeline for the H. Sapiens proteome subset.

Usage
-----
    python src/process_proteome.py [--n N] [--seed SEED] [--t-p T_P]
                                   [--pdb-dir DIR] [--out PATH]

Workflow
--------
1. Fetch all reviewed H. Sapiens UniProt accessions from the UniProt REST API
   (cached to data/results/human_accessions.txt on first run).
2. Sample N accessions at random (default 500, seed 42).
3. Download the AlphaFold-DB PDB for each (current latest version, ~v6).
   Downloads run in parallel (10 workers); failures are logged and skipped.
4. For each downloaded PDB: parse pLDDT, run build_pd_and_ncc, compute
   f⁺_cp, mean_persistence, H_p, ncc_max, PLM(t_p).
5. Write a tidy Parquet file to --out (default data/results/proteome_subset_N500.parquet).

Each AlphaFold-DB fragment file (F1, F2, …) is treated as a separate structure.
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

# Ensure src/ is on the path when run directly.
_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from download import download_structure
from filtration import build_pd_and_ncc
from parse import parse_alphafold_pdb
from plm import plm as compute_plm
from statistics import f_cp_plus, mean_persistence, ncc_max, persistence_entropy

log = logging.getLogger(__name__)

_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_PDB_DIR = _ROOT / "data" / "hsapiens"
_DEFAULT_OUT = _ROOT / "data" / "results" / "proteome_subset_N500.parquet"
_ACCESSION_CACHE = _ROOT / "data" / "results" / "human_accessions.txt"

_UNIPROT_STREAM = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id%3A9606+AND+reviewed%3Atrue"
    "&format=list"
)


# ---------------------------------------------------------------------------
# Step 1: accession list
# ---------------------------------------------------------------------------

def fetch_human_accessions(cache_path: Path = _ACCESSION_CACHE) -> list[str]:
    """
    Return all reviewed H. Sapiens UniProt accessions.

    Downloads from UniProt on first call and caches to disk; subsequent
    calls read from the cache.
    """
    if cache_path.exists():
        log.info("Loading accessions from cache: %s", cache_path)
        return cache_path.read_text().splitlines()

    log.info("Fetching H. Sapiens accessions from UniProt…")
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(_UNIPROT_STREAM, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        accessions = [line.strip() for line in resp.iter_lines(decode_unicode=True)
                      if line.strip()]

    cache_path.write_text("\n".join(accessions))
    log.info("Fetched %d accessions → %s", len(accessions), cache_path)
    return accessions


# ---------------------------------------------------------------------------
# Step 2: download subset
# ---------------------------------------------------------------------------

def _download_one(uid: str, pdb_dir: Path) -> tuple[str, Path | None]:
    """Download a single PDB; return (uid, path) or (uid, None) on failure."""
    try:
        path = download_structure(uid, pdb_dir)
        return uid, path
    except Exception as exc:
        log.warning("Skipping %s: %s", uid, exc)
        return uid, None


def download_subset(
    accessions: list[str],
    n: int,
    seed: int,
    pdb_dir: Path,
    workers: int = 10,
) -> list[Path]:
    """
    Sample *n* accessions at random and download their AlphaFold-DB PDB files.

    Returns the list of successfully downloaded paths.
    """
    rng = np.random.default_rng(seed)
    sample = rng.choice(accessions, size=min(n, len(accessions)), replace=False).tolist()

    pdb_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_download_one, uid, pdb_dir): uid for uid in sample}
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="Downloading", unit="pdb"):
            uid, path = future.result()
            if path is not None:
                paths.append(path)

    log.info("Downloaded %d / %d structures.", len(paths), len(sample))
    return paths


# ---------------------------------------------------------------------------
# Step 3: per-protein statistics
# ---------------------------------------------------------------------------

def _compute_stats(pdb_path: Path, t_p: float) -> dict | None:
    """Parse a single PDB and compute all five statistics.  Returns None on error."""
    try:
        data = parse_alphafold_pdb(pdb_path)
    except Exception as exc:
        log.warning("Parse error %s: %s", pdb_path.name, exc)
        return None

    plddt = data["plddt"]
    n = len(plddt)
    if n < 2:
        return None

    try:
        pd_arr, ncc_curve = build_pd_and_ncc(plddt)
        return {
            "uniprot_id":       data["uniprot_id"],
            "filename":         pdb_path.name,
            "n_residues":       n,
            "f_cp_plus":        f_cp_plus(pd_arr),
            "mean_persistence": mean_persistence(pd_arr),
            "H_p":              persistence_entropy(pd_arr),
            "ncc_max":          ncc_max(ncc_curve),
            "PLM":              compute_plm(ncc_curve, n_residues=n, t_p=t_p),
        }
    except Exception as exc:
        log.warning("Filtration error %s: %s", pdb_path.name, exc)
        return None


def run_pipeline(
    pdb_dir: Path,
    out_path: Path,
    t_p: float = 0.025,
) -> pd.DataFrame:
    """
    Process all PDB files in *pdb_dir*, write results to *out_path* (Parquet).

    Returns the resulting DataFrame.
    """
    pdbs = sorted(pdb_dir.glob("AF-*-model_v*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No AlphaFold PDB files found in {pdb_dir}")

    rows = []
    for pdb in tqdm(pdbs, desc="Computing stats", unit="pdb"):
        row = _compute_stats(pdb, t_p)
        if row is not None:
            rows.append(row)

    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    log.info("Wrote %d rows → %s", len(df), out_path)
    return df


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(description="H. Sapiens proteome subset pipeline")
    parser.add_argument("--n",       type=int,   default=500,
                        help="Number of proteins to sample (default 500)")
    parser.add_argument("--seed",    type=int,   default=42,
                        help="Random seed for sampling (default 42)")
    parser.add_argument("--t-p",    type=float, default=0.025,
                        help="PLM relative persistence threshold (default 0.025)")
    parser.add_argument("--pdb-dir", type=Path,  default=_DEFAULT_PDB_DIR,
                        help="Directory for downloaded PDB files")
    parser.add_argument("--out",     type=Path,  default=None,
                        help="Output Parquet path (default auto-named by N)")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip download; process PDBs already in --pdb-dir")
    args = parser.parse_args()

    out_path = args.out or (_ROOT / "data" / "results" / f"proteome_subset_N{args.n}.parquet")

    if not args.skip_download:
        accessions = fetch_human_accessions()
        print(f"Total H. Sapiens accessions: {len(accessions)}")
        download_subset(accessions, n=args.n, seed=args.seed, pdb_dir=args.pdb_dir)

    df = run_pipeline(args.pdb_dir, out_path, t_p=args.t_p)
    print(f"\nProcessed {len(df)} proteins.")
    print(df[["n_residues", "f_cp_plus", "H_p", "PLM"]].describe().round(3).to_string())

    from scipy.stats import pearsonr
    r, pval = pearsonr(df["f_cp_plus"], df["H_p"])
    print(f"\nPearson r(f⁺_cp, H_p) = {r:.4f}  (p = {pval:.2e})")
    print(f"Target: r ≈ 0.97  ({'✓ PASS' if r >= 0.85 else '✗ below threshold'})")
