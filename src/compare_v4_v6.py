"""
compare_v4_v6.py — compare proteome-scale results between AlphaFold-DB v4 and v6.

Usage
-----
    python src/compare_v4_v6.py

Reads:
    data/results/proteome_full.parquet      (v4)
    data/results/proteome_v6_full.parquet   (v6)

Reports:
    - Protein counts
    - Median / mean f⁺_cp, H_p, PLM
    - Pearson r(f⁺_cp, H_p) for each version
    - Fragmentation candidates: n≥200, H_p≥0.25, PLM≥3 at t_p=0.025
    - Per-column shift between v4 and v6
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

_ROOT  = Path(__file__).resolve().parent.parent
_V4    = _ROOT / "data" / "results" / "proteome_full.parquet"
_V6    = _ROOT / "data" / "results" / "proteome_v6_full.parquet"


def load(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    # normalise column name capitalisation across runs
    df.columns = [c.lower() for c in df.columns]
    if "h_p" not in df.columns and "hp" in df.columns:
        df = df.rename(columns={"hp": "h_p"})
    return df


def report(label: str, df: pd.DataFrame) -> None:
    print(f"\n{'='*60}")
    print(f"  {label}  (N = {len(df):,})")
    print(f"{'='*60}")

    for col in ("f_cp_plus", "h_p", "ncc_max", "plm"):
        if col not in df.columns:
            continue
        s = df[col]
        print(f"  {col:18s}  median={s.median():.4f}  mean={s.mean():.4f}"
              f"  std={s.std():.4f}")

    r, pval = pearsonr(df["f_cp_plus"], df["h_p"])
    print(f"\n  Pearson r(f⁺_cp, H_p) = {r:.4f}  (p = {pval:.2e})")
    print(f"  Paper target           ≈ 0.97")

    cands = df.query("n_residues >= 200 and h_p >= 0.25 and plm >= 3")
    print(f"\n  Fragmentation candidates (n≥200, H_p≥0.25, PLM≥3): {len(cands):,}")
    print(f"  Paper target: ≈ 86")


def compare(v4: pd.DataFrame, v6: pd.DataFrame) -> None:
    print(f"\n{'='*60}")
    print("  v4 → v6  shift (median values)")
    print(f"{'='*60}")
    for col in ("f_cp_plus", "h_p", "ncc_max"):
        if col not in v4.columns or col not in v6.columns:
            continue
        delta = v6[col].median() - v4[col].median()
        print(f"  Δ {col:18s}  {delta:+.4f}")

    # protein count difference
    in_v6_not_v4 = set(v6["uniprot_id"]) - set(v4["uniprot_id"])
    in_v4_not_v6 = set(v4["uniprot_id"]) - set(v6["uniprot_id"])
    print(f"\n  Proteins in v6 only : {len(in_v6_not_v4):,}")
    print(f"  Proteins in v4 only : {len(in_v4_not_v6):,}")


if __name__ == "__main__":
    if not _V4.exists():
        raise FileNotFoundError(f"v4 parquet not found: {_V4}")
    if not _V6.exists():
        raise FileNotFoundError(f"v6 parquet not found: {_V6}")

    v4 = load(_V4)
    v6 = load(_V6)

    report("AlphaFold-DB v4", v4)
    report("AlphaFold-DB v6", v6)
    compare(v4, v6)
