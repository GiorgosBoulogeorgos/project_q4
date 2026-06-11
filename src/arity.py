"""
arity.py — Cα packing-arity descriptors for Q1 cross-correlation.

Arity of residue i: number of *other* Cα atoms within 10 Å.
Arity signature of a protein: (arity_25, arity_75) — the 25th and 75th
percentiles of the per-residue arity distribution.

Reference: Cazals & Sarti (2025), Q1 — 3D packing analysis and arity maps.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def arity_per_residue(ca_coords: np.ndarray, r: float = 10.0) -> np.ndarray:
    """
    Count Cα neighbours within distance r Å for each residue (self excluded).

    Parameters
    ----------
    ca_coords : (n, 3) float array of Cα coordinates in Angstroms.
    r         : neighbourhood radius in Å (default 10.0, paper convention).

    Returns
    -------
    (n,) int array — arity of each residue.
    """
    if len(ca_coords) == 0:
        return np.array([], dtype=np.int32)
    tree = cKDTree(ca_coords)
    counts = tree.query_ball_point(ca_coords, r=r, return_length=True)
    return np.asarray(counts, dtype=np.int32) - 1  # subtract self


def arity_signature(
    ca_coords: np.ndarray,
    r: float = 10.0,
    q1: float = 0.25,
    q2: float = 0.75,
) -> tuple[int, int]:
    """
    Return the (arity_q1, arity_q2) signature of a protein.

    Default quantiles: 25th and 75th (paper convention).

    Per Cazals & Sarti (2025), Def. 2 / Eq. 2, the arity at quantile q is the
    *inverse empirical CDF*: the smallest observed arity a such that
    CDF_Cα(a) >= q.  This is ``np.quantile(..., method="inverted_cdf")`` — it
    always returns an arity that actually occurs in the protein.  (Plain
    ``np.percentile`` uses linear interpolation, which can return a value
    between two observed arities — e.g. 9 for a distribution whose arities jump
    3 -> 12 — and therefore does not match the paper's definition.)
    """
    arities = arity_per_residue(ca_coords, r=r)
    return (
        int(np.quantile(arities, q1, method="inverted_cdf")),
        int(np.quantile(arities, q2, method="inverted_cdf")),
    )
