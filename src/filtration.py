"""
filtration.py — path-graph pLDDT filtration with Union-Find (Elder Rule).

A polypeptide of n residues is modelled as a path graph
0 — 1 — 2 — … — (n−1); each vertex carries a pLDDT score.  The filtration
inserts residues in DECREASING pLDDT order; edges are activated when both
endpoints are present.  This is the standard sublevel-set filtration described
in Cazals & Sarti (2025), using u = −pLDDT as the filter function.

Canonical function
------------------
``build_pd_and_ncc(plddt)``
    One residue per insertion step, raw pLDDT values (no discretisation).
    Residues are sorted by decreasing pLDDT; ties are broken by ascending
    index (stable sort).  This is the paper-faithful implementation.

Elder Rule
----------
When two components merge, the one with HIGHER birth_pLDDT is "older" and
survives as root; the younger one dies with death = current_pLDDT.
Tie-break: lower root index survives.

Returns
-------
pd : np.ndarray, shape (n−1, 2)
    Persistence diagram — each row is (birth_pLDDT, death_pLDDT) with
    birth_pLDDT >= death_pLDDT.  Accretions have birth == death.
ncc_curve : np.ndarray, shape (n, 2)
    One row per insertion step: (pLDDT_value, N_cc_after_this_insertion).
    Ordered by decreasing pLDDT (insertion order).
"""

from __future__ import annotations

import numpy as np


def build_pd_and_ncc(
    plddt: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the path-graph pLDDT filtration (paper-faithful implementation).

    Inserts one residue per step in decreasing pLDDT order.  Edges are
    activated as soon as both endpoints are present.  Ties in pLDDT are
    broken by ascending residue index (stable sort).

    Parameters
    ----------
    plddt:
        1-D array of per-residue pLDDT values (0–100 scale), in residue order.

    Returns
    -------
    pd:
        Array of shape (n−1, 2) with columns [birth_pLDDT, death_pLDDT].
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT_value, N_cc].
        One row per insertion step, in decreasing-pLDDT order.
    """
    plddt = np.asarray(plddt, dtype=np.float64)
    n = len(plddt)
    if n == 0:
        return np.empty((0, 2), dtype=np.float64), np.empty((0, 2), dtype=np.float64)

    # Stable sort preserves left-to-right order among ties.
    order = np.argsort(-plddt, kind="stable")

    parent = np.arange(n, dtype=np.intp)
    birth_of_root = np.zeros(n, dtype=np.float64)
    inserted = np.zeros(n, dtype=bool)

    pd_pairs: list[tuple[float, float]] = []
    ncc_curve: list[tuple[float, int]] = []
    n_cc = 0

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]  # path halving
            i = parent[i]
        return int(i)

    def _merge(a: int, b: int, death: float) -> None:
        nonlocal n_cc
        ra, rb = find(a), find(b)
        if ra == rb:
            pd_pairs.append((death, death))
            return
        ba, bb = birth_of_root[ra], birth_of_root[rb]
        if ba > bb or (ba == bb and ra < rb):
            parent[rb] = ra
            pd_pairs.append((bb, death))
        else:
            parent[ra] = rb
            pd_pairs.append((ba, death))
        n_cc -= 1

    for idx in order:
        i = int(idx)
        p = float(plddt[i])

        parent[i] = i
        birth_of_root[i] = p
        inserted[i] = True
        n_cc += 1

        if i > 0 and inserted[i - 1]:
            _merge(i, i - 1, p)
        if i < n - 1 and inserted[i + 1]:
            _merge(i, i + 1, p)

        ncc_curve.append((p, n_cc))

    pd = np.array(pd_pairs, dtype=np.float64).reshape(-1, 2)
    return pd, np.array(ncc_curve, dtype=np.float64)
