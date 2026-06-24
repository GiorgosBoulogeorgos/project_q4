"""
compare_organisms.py — cross-organism comparison of the proteome-scale Q4
statistics (AlphaFold-DB v6, integer-pLDDT batched pipeline).

Usage
-----
    python src/compare_organisms.py

Reads (all produced by src/run_organism.py / the v6 H. sapiens run):
    data/results/proteome_v6_full.parquet            (H. sapiens)
    data/results/proteome_<org>_v6_full.parquet      (mouse / rat / yeast)
    data/results/proteome_<org>_v6_arity.parquet     (arity signatures)

Reports a single apples-to-apples table — N, Pearson r(f⁺_cp, H_p), median
f⁺_cp / H_p, max H_p, Figure-8 fragmentation candidates, and the arity
enrichment of those candidates. The H. sapiens baseline is computed *here*
from our own run, not taken from the paper, so the comparison is on the
identical pipeline. The paper's H. sapiens targets (r ≈ 0.97, ≈ 86
candidates) are printed only as an external reference.

Also writes a bar chart of Pearson r per organism to
    data/results/cross_organism_pearson.png / .pdf
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, pearsonr

_ROOT    = Path(__file__).resolve().parent.parent
_RESULTS = _ROOT / "data" / "results"          # input parquets
_FIG_DIR = _ROOT / "results" / "figures"       # report figures (report/figures -> here)

# label → (full parquet, arity parquet). All v6, so arity merges cleanly on
# (uniprot_id, fragment) within each organism.
_ORGANISMS: dict[str, tuple[Path, Path]] = {
    "H. sapiens":    (_RESULTS / "proteome_v6_full.parquet",
                      _RESULTS / "proteome_v6_arity.parquet"),
    "M. musculus":   (_RESULTS / "proteome_mouse_v6_full.parquet",
                      _RESULTS / "proteome_mouse_v6_arity.parquet"),
    "R. norvegicus": (_RESULTS / "proteome_rat_v6_full.parquet",
                      _RESULTS / "proteome_rat_v6_arity.parquet"),
    "D. melanogaster": (_RESULTS / "proteome_drome_v6_full.parquet",
                        _RESULTS / "proteome_drome_v6_arity.parquet"),
    "S. cerevisiae": (_RESULTS / "proteome_yeast_v6_full.parquet",
                      _RESULTS / "proteome_yeast_v6_arity.parquet"),
    "P. aeruginosa": (_RESULTS / "proteome_pseae_v6_full.parquet",
                      _RESULTS / "proteome_pseae_v6_arity.parquet"),
}

_PAPER_R    = 0.97   # Cazals & Sarti, H. sapiens (Fig. 7)
_PAPER_CAND = 86     # Cazals & Sarti, H. sapiens (Fig. 8)


def load(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df.columns = [c.lower() for c in df.columns]
    return df


def arity_enrichment(full: pd.DataFrame, arity: pd.DataFrame):
    """
    Replicate the notebook's Q1×Q4 cross-correlation: find the median arity
    bin of the fragmentation candidates and test its enrichment against the
    whole proteome (one-sided Fisher exact).

    Returns (enrichment, p_value, (a25, a75), n_candidates) or None if the
    merge yields nothing.
    """
    merged = full.merge(
        arity[["uniprot_id", "fragment", "arity_25", "arity_75"]],
        on=["uniprot_id", "fragment"],
    )
    if merged.empty:
        return None
    cands = merged.query("n_residues >= 200 and h_p >= 0.25 and plm >= 3")
    if len(cands) == 0:
        return None

    a25 = int(round(cands["arity_25"].median()))
    a75 = int(round(cands["arity_75"].median()))
    in_bin_c  = int(((cands["arity_25"]  == a25) & (cands["arity_75"]  == a75)).sum())
    in_bin_bg = int(((merged["arity_25"] == a25) & (merged["arity_75"] == a75)).sum())
    n_c, n_bg = len(cands), len(merged)
    table = [[in_bin_c, n_c - in_bin_c],
             [in_bin_bg - in_bin_c, n_bg - n_c - (in_bin_bg - in_bin_c)]]
    _, pval = fisher_exact(table, alternative="greater")
    enrichment = (in_bin_c / n_c) / (in_bin_bg / n_bg) if in_bin_bg else float("nan")
    return enrichment, pval, (a25, a75), n_c


def summarise(label: str, full_path: Path, arity_path: Path) -> dict | None:
    if not full_path.exists():
        print(f"  [skip] {label}: missing {full_path.name}")
        return None
    df = load(full_path)
    r, _ = pearsonr(df["f_cp_plus"], df["h_p"])
    cands = df.query("n_residues >= 200 and h_p >= 0.25 and plm >= 3")

    row = {
        "organism":  label,
        "n":         len(df),
        "r":         r,
        "med_fcp":   df["f_cp_plus"].median(),
        "med_hp":    df["h_p"].median(),
        "max_hp":    df["h_p"].max(),
        "n_cand":    len(cands),
        "enrich":    None,
        "enrich_p":  None,
        "bin":       None,
    }
    if arity_path.exists():
        res = arity_enrichment(df, load(arity_path))
        if res is not None:
            row["enrich"], row["enrich_p"], row["bin"], _ = res
    return row


def print_table(rows: list[dict]) -> None:
    # cells joined with a 2-space gap so columns never collide
    sep = "  "
    fmt = lambda cells: sep.join(cells)
    head = fmt([f"{'Organism':<14}", f"{'N':>8}", f"{'r':>7}",
                f"{'med f+_cp':>9}", f"{'med H_p':>7}", f"{'max H_p':>7}",
                f"{'Fig8':>5}", f"{'arity enrich':<18}"])
    rule = "=" * len(head)
    print("\n" + rule)
    print("  Cross-organism comparison (AlphaFold-DB v6, integer+batched)")
    print(rule)
    print(head)
    print("-" * len(head))
    for x in rows:
        enr = (f"{x['enrich']:.1f}x (p={x['enrich_p']:.1e})"
               if x["enrich"] is not None else "n/a")
        print(fmt([f"{x['organism']:<14}", f"{x['n']:>8,}", f"{x['r']:>7.4f}",
                   f"{x['med_fcp']:>9.4f}", f"{x['med_hp']:>7.4f}",
                   f"{x['max_hp']:>7.4f}", f"{x['n_cand']:>5d}", f"{enr:<18}"]))
    print("-" * len(head))
    print(f"  External reference (paper, H. sapiens): r ≈ {_PAPER_R}, "
          f"Fig8 ≈ {_PAPER_CAND} candidates")


def plot_pearson(rows: list[dict], out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                       # pragma: no cover
        print(f"  [plot skipped] matplotlib unavailable: {exc}")
        return

    labels = [r["organism"] for r in rows]
    vals   = [r["r"] for r in rows]
    human  = next((r["r"] for r in rows if r["organism"] == "H. sapiens"), None)

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    ax.bar(labels, vals, color="#4C72B0")
    if human is not None:
        ax.axhline(human, ls="--", c="#444",
                   label=f"our H. sapiens = {human:.3f}")
    ax.axhline(_PAPER_R, ls=":", c="#C44E52", label=f"paper ≈ {_PAPER_R}")
    ax.set_ylim(0.7, 1.0)
    ax.set_ylabel(r"Pearson $r(f^+_{cp}, H_p)$")
    ax.set_title("Cross-organism correlation (AlphaFold-DB v6)")
    ax.legend(fontsize=8, loc="upper right")
    # italicise species names and rotate so 6 labels don't collide
    # (build mathtext outside the f-string: py3.9 forbids backslashes in f-strings)
    italic = [r"$\it{" + l.replace(" ", r"\ ") + "}$" for l in labels]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(italic, rotation=20, ha="right", fontsize=8)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)
    print(f"\n  Wrote {out.with_suffix('.png')} (+ .pdf)")


if __name__ == "__main__":
    rows = [r for r in (summarise(lbl, f, a)
                        for lbl, (f, a) in _ORGANISMS.items()) if r]
    if not rows:
        raise SystemExit("No organism parquets found under data/results/.")
    print_table(rows)
    plot_pearson(rows, _FIG_DIR / "cross_organism_pearson")
