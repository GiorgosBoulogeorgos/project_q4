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

    Default percentiles: 25th and 75th (paper convention).
    Values are floored to int via np.percentile's linear interpolation.
    """
    arities = arity_per_residue(ca_coords, r=r)
    return (
        int(np.percentile(arities, q1 * 100)),
        int(np.percentile(arities, q2 * 100)),
    )
