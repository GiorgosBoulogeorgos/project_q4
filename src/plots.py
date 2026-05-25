"""
plots.py — Figure generation for Q4 proteome analysis.

Functions
---------
figure7(df, out_path, ...)
    Scatter of f⁺_cp vs H_p (reproduces Figure 7 of the paper).
figure8_scatter(df, out_path, ...)
    3D scatter of n_residues × H_p × f⁺_cp (reproduces Figure 8A).
render_structure(pdb_path, out_path, title='')
    Matplotlib 3D Cα trace coloured by pLDDT (approximates Fig. 8 B/C/D).
figure8_panels(fragmented_df, out_dir, n_examples=4)
    Download and render the top-PLM fragmented proteins.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


# ---------------------------------------------------------------------------
# Prototype constants used for overlay on all figures
# ---------------------------------------------------------------------------

_PROTOTYPES = [
    {"uid": "P15121",     "label": "P15121\n(Ordered)",     "color": "#1f77b4", "marker": "*"},
    {"uid": "A0A0G2L439", "label": "A0A0G2L439\n(Disordered)", "color": "#d62728", "marker": "*"},
    {"uid": "Q9VQS4",     "label": "Q9VQS4\n(Mixed)",       "color": "#2ca02c", "marker": "*"},
]


# ---------------------------------------------------------------------------
# Figure 7 — f⁺_cp vs H_p scatter
# ---------------------------------------------------------------------------

def figure7(
    df: pd.DataFrame,
    out_path: Path,
    proto_rows: list[dict] | None = None,
    title: str = "",
) -> None:
    """
    Reproduce Figure 7: scatter of f⁺_cp vs H_p, coloured by log(n_residues).

    Parameters
    ----------
    df:
        DataFrame with columns f_cp_plus, H_p, n_residues, uniprot_id.
    out_path:
        Where to save the PNG.
    proto_rows:
        List of dicts with keys uid, f_cp_plus, H_p, label, color.
        If None, no prototype overlay is drawn.
    title:
        Optional title prefix.
    """
    r, _ = pearsonr(df["f_cp_plus"], df["H_p"])

    fig, ax = plt.subplots(figsize=(8, 7))

    # Main scatter — colour by log10(n_residues)
    n_res = df["n_residues"].values.clip(1)
    sc = ax.scatter(
        df["f_cp_plus"], df["H_p"],
        c=np.log10(n_res), cmap="viridis",
        s=12, alpha=0.45, linewidths=0,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("log₁₀(n residues)", fontsize=10)
    cbar.set_ticks([2, 2.5, 3, 3.5])
    cbar.set_ticklabels(["100", "316", "1 000", "3 162"])

    # Linear fit
    xg = np.linspace(df["f_cp_plus"].min(), df["f_cp_plus"].max(), 200)
    m, b = np.polyfit(df["f_cp_plus"], df["H_p"], 1)
    ax.plot(xg, m * xg + b, "k--", lw=1.2, alpha=0.7,
            label=f"r = {r:.4f}")

    # Prototype overlays
    if proto_rows:
        for p in proto_rows:
            ax.scatter(p["f_cp_plus"], p["H_p"], s=200,
                       marker="*", color=p["color"], zorder=6,
                       edgecolors="k", linewidths=0.5,
                       label=p["label"].replace("\n", " "))
            ax.annotate(
                p["uid"],
                (p["f_cp_plus"], p["H_p"]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=7.5, color=p["color"], fontweight="bold",
            )

    # Annotate extremes
    for side in ["max", "min"]:
        idx = df["H_p"].idxmax() if side == "max" else df["H_p"].idxmin()
        row = df.loc[idx]
        ax.annotate(
            f"{row['uniprot_id']} ({side} H_p)",
            (row["f_cp_plus"], row["H_p"]),
            xytext=(8, -10), textcoords="offset points",
            fontsize=6.5, color="gray",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.6),
        )

    ax.set_xlabel("f⁺_cp", fontsize=13)
    ax.set_ylabel("H_p",   fontsize=13)
    n_prot = len(df)
    full_title = f"H. Sapiens — {n_prot:,} AlphaFold predictions\n"
    full_title += f"Pearson r(f⁺_cp, H_p) = {r:.4f}  [paper target ≈ 0.97]"
    if title:
        full_title = title + "\n" + full_title
    ax.set_title(full_title, fontsize=11)
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.85)
    ax.grid(alpha=0.2)

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
    print(f"Saved Figure 7 → {out_path}")


# ---------------------------------------------------------------------------
# Figure 8A — 3D scatter
# ---------------------------------------------------------------------------

def figure8_scatter(
    df: pd.DataFrame,
    out_path: Path,
    strict_mask: pd.Series | None = None,
    calibrated_mask: pd.Series | None = None,
) -> None:
    """
    Reproduce Figure 8A: 3D scatter of n_residues × H_p × f⁺_cp.

    Fragmented proteins (calibrated filter) are highlighted in red.
    Strict-threshold proteins (paper's H_p≥0.25) are highlighted in orange.
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    x = df["n_residues"].values
    y = df["H_p"].values
    z = df["f_cp_plus"].values

    bg_mask = np.ones(len(df), dtype=bool)
    if calibrated_mask is not None:
        bg_mask &= ~calibrated_mask.values

    ax.scatter(x[bg_mask], y[bg_mask], z[bg_mask],
               c="steelblue", s=4, alpha=0.20, linewidths=0,
               label=f"background (n={bg_mask.sum():,})")

    if calibrated_mask is not None and calibrated_mask.any():
        n_frag = calibrated_mask.sum()
        ax.scatter(x[calibrated_mask], y[calibrated_mask], z[calibrated_mask],
                   c="crimson", s=30, alpha=0.85, linewidths=0,
                   label=f"fragmented H_p≥0.40, PLM≥3 (n={n_frag})")

    if strict_mask is not None and strict_mask.any():
        n_strict = strict_mask.sum()
        ax.scatter(x[strict_mask], y[strict_mask], z[strict_mask],
                   c="orange", s=40, alpha=0.9, marker="D", linewidths=0,
                   label=f"paper H_p≥0.25, PLM≥3 (n={n_strict})")

    ax.set_xlabel("n residues", fontsize=10, labelpad=8)
    ax.set_ylabel("H_p",        fontsize=10, labelpad=8)
    ax.set_zlabel("f⁺_cp",      fontsize=10, labelpad=8)
    ax.set_title("Figure 8A — Fragmented AlphaFold predictions\n"
                 "H. Sapiens  ·  filter: n≥200, PLM≥3", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.view_init(elev=20, azim=-60)

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
    print(f"Saved Figure 8A → {out_path}")


# ---------------------------------------------------------------------------
# Structure rendering — Cα trace coloured by pLDDT
# ---------------------------------------------------------------------------

def render_structure(
    pdb_path: Path,
    out_path: Path,
    title: str = "",
) -> None:
    """
    Render a Cα trace with pLDDT colouring and save to PNG.

    Colour scheme mirrors AlphaFold-DB:
        [0, 50)    orange-red   (very low confidence)
        [50, 70)   yellow       (low confidence)
        [70, 90)   cyan         (confident)
        [90, 100]  dark-blue    (very high confidence)
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from parse import parse_alphafold_pdb

    data   = parse_alphafold_pdb(pdb_path)
    coords = data["ca_coords"]
    plddt  = data["plddt"]
    uid    = data["uniprot_id"]
    n      = len(plddt)

    fig = plt.figure(figsize=(6, 5))
    ax  = fig.add_subplot(111, projection="3d")

    # Backbone line
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
            color="lightgray", lw=0.6, alpha=0.5, zorder=1)

    # Cα atoms coloured by pLDDT
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                    c=plddt, cmap="RdYlBu_r", vmin=0, vmax=100,
                    s=10, alpha=0.85, linewidths=0, zorder=2)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.05)
    cbar.set_label("pLDDT", fontsize=9)
    cbar.set_ticks([0, 50, 70, 90, 100])

    ax.set_xlabel("x (Å)", fontsize=7)
    ax.set_ylabel("y (Å)", fontsize=7)
    ax.set_zlabel("z (Å)", fontsize=7)
    ax.tick_params(labelsize=6)

    head = title or f"{uid}  (n={n}  mean pLDDT={plddt.mean():.1f})"
    ax.set_title(head, fontsize=9)
    ax.grid(alpha=0.15)

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 8 B/C/D panels — download and render top-PLM fragmented proteins
# ---------------------------------------------------------------------------

def figure8_panels(
    fragmented_df: pd.DataFrame,
    out_dir: Path,
    n_examples: int = 4,
) -> list[Path]:
    """
    Download and render Cα traces for the top-PLM fragmented proteins.

    Returns the list of saved PNG paths.
    """
    from download import download_structure

    top = (fragmented_df
           .sort_values(["PLM", "H_p"], ascending=False)
           .drop_duplicates("uniprot_id")
           .head(n_examples))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    for _, row in top.iterrows():
        uid = row["uniprot_id"]
        try:
            pdb_path = download_structure(uid, out_dir)
        except Exception as exc:
            print(f"  Skipping {uid}: {exc}")
            continue
        out_png = out_dir / f"{uid}_structure.png"
        title = (f"{uid}  n={int(row['n_residues'])}  "
                 f"H_p={row['H_p']:.3f}  PLM={int(row['PLM'])}")
        render_structure(pdb_path, out_png, title=title)
        saved.append(out_png)
        print(f"  Rendered {uid} → {out_png.name}")

    return saved
