"""
filtration.py — path-graph pLDDT filtration with Union-Find (Elder Rule).

A polypeptide of n residues is modelled as a path graph
0 — 1 — 2 — … — (n−1); each vertex carries a pLDDT score.  The filtration
inserts residues in DECREASING pLDDT order; edges are activated when both
endpoints are present.

Canonical function
------------------
``build_pd_and_ncc(plddt, discretize=np.round)``
    Batch insertion: all residues that share the same *discretised* pLDDT
    level are inserted simultaneously.  This is the production function used
    throughout the pipeline.

    Motivation: AlphaFold-DB v4 stores high-precision float pLDDT values
    (up to 108 distinct values for a 316-residue protein).  The paper's
    analysis appears to have used effectively integer-valued pLDDT, which
    collapses many residues to the same level.  Batch insertion handles ties
    principally: all merges at level L produce PD pairs with birth = death = L
    (accretions), while cross-level merges produce positive-persistence pairs.
    The result is a much smaller |P|, consistent with the paper's H_p targets.

    The ``discretize`` keyword (default ``np.round``) lets callers override
    the rounding strategy without touching any production code.  Pass
    ``np.floor`` for an alternative discretisation (Experiment C in the
    diagnostic run).

Reference / ablation function
------------------------------
``build_pd_and_ncc_single(plddt)``
    One residue inserted per step, using the raw (undiscretised) pLDDT
    values.  Kept for reference and ablation.  Produces H_p values
    0.14–0.33 above the paper's Figure 3 targets for the three prototype
    proteins; the discrepancy is attributed to the higher float precision of
    v4/v6 AlphaFold-DB files.

    (Renamed from the original ``build_pd_and_ncc``.  The canonical name now
    belongs to the batch variant.)

Size-elder ablation
-------------------
``build_pd_and_ncc_batch_size_elder(plddt, discretize=np.round)``
    Same as the canonical function except when two components share the same
    birth_pLDDT: the larger component survives instead of the lower-index one.
    Produces identical PD pairs to the canonical function (within-batch merges
    always have birth = death = L regardless of which component survives).
    Kept for completeness.

Elder Rule (common to all variants)
------------------------------------
When two components merge, the one with HIGHER birth_pLDDT is "older" and
survives as root; the younger one dies with death = current_pLDDT.
Tie-break: lower root index survives (canonical) or larger size (size-elder).

Returns (all variants)
----------------------
pd : np.ndarray, shape (n−1, 2)
    Persistence diagram — each row is (birth_pLDDT, death_pLDDT) with
    birth_pLDDT >= death_pLDDT.  Accretions have birth == death.
ncc_curve : np.ndarray, shape (n, 2)
    One row per residue: (pLDDT_level, N_cc_after_batch).
    For the batch variant, all residues in a batch share the same N_cc value
    (post-merge N_cc for that level).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Canonical function (batch insertion, discretised pLDDT)
# ---------------------------------------------------------------------------

def build_pd_and_ncc(
    plddt: np.ndarray,
    discretize=np.round,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the path-graph pLDDT filtration with batch insertion.

    All residues that share the same discretised pLDDT level are inserted
    simultaneously.  Edges within a batch produce accretion pairs
    (birth = death = level); edges crossing two different levels produce
    positive-persistence pairs.

    Parameters
    ----------
    plddt:
        1-D array of per-residue pLDDT values (0–100 scale), in residue order.
    discretize:
        Callable applied to the pLDDT array before the filtration.
        Default: ``np.round`` (integer levels).
        Alternative: ``np.floor`` for floor-discretisation.

    Returns
    -------
    pd:
        Array of shape (n−1, 2) with columns [birth_pLDDT, death_pLDDT].
        Birth and death are the discretised integer levels cast to float.
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT_level, N_cc].
        All residues in the same batch share the same post-merge N_cc value.
    """
    plddt = np.asarray(plddt, dtype=np.float64)
    n = len(plddt)
    if n == 0:
        return np.empty((0, 2), dtype=np.float64), np.empty((0, 2), dtype=np.float64)

    plddt_disc = discretize(plddt).astype(np.int64)
    levels = np.unique(plddt_disc)[::-1]  # process highest level first

    parent = np.arange(n, dtype=np.intp)
    birth_of_root = np.zeros(n, dtype=np.float64)
    inserted = np.zeros(n, dtype=bool)

    pd_pairs: list[tuple[float, float]] = []
    ncc_rows: list[tuple[float, int]] = []
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

    for level in levels:
        L = float(level)
        batch = np.where(plddt_disc == level)[0]  # ascending index order

        # Step 1: register all residues in this batch as singletons.
        for idx in batch:
            i = int(idx)
            parent[i] = i
            birth_of_root[i] = L
            inserted[i] = True
            n_cc += 1

        # Step 2: activate edges — each path-graph edge processed exactly once.
        # Left-check (always): covers within-batch edges (right residue checks
        #   its left neighbour, which was already inserted in this batch).
        # Right-check (only for cross-level right neighbours): covers edges
        #   from this batch to already-inserted higher-level residues on the
        #   right; within-batch right edges are deliberately skipped here and
        #   handled by the left-check of the right residue.
        for idx in batch:
            i = int(idx)
            if i > 0 and inserted[i - 1]:
                _merge(i, i - 1, L)
            if i < n - 1 and inserted[i + 1] and plddt_disc[i + 1] > level:
                _merge(i, i + 1, L)

        # One ncc_curve row per residue; all share the post-merge N_cc.
        for _ in batch:
            ncc_rows.append((L, n_cc))

    pd_arr = np.array(pd_pairs, dtype=np.float64).reshape(-1, 2)
    return pd_arr, np.array(ncc_rows, dtype=np.float64)


# ---------------------------------------------------------------------------
# Reference / ablation function (single insertion, raw float pLDDT)
# ---------------------------------------------------------------------------

def build_pd_and_ncc_single(plddt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the path-graph pLDDT filtration, inserting one residue per step.

    Uses raw (undiscretised) pLDDT values.  Kept for reference and ablation;
    produces H_p values 0.14–0.33 above the paper's Figure 3 targets for
    current AlphaFold-DB data.

    Three event types at each insertion step:
      A) 0 present neighbours  → N_cc += 1;  no PD pair.
      B) 1 present neighbour   → N_cc unchanged; accretion PD pair (b, b).
      C) 2 present neighbours  → N_cc −= 1;  1 accretion + 1 positive pair.

    Parameters
    ----------
    plddt:
        1-D array of per-residue pLDDT values (0–100 scale), in residue order.

    Returns
    -------
    pd:
        Array of shape (n−1, 2) with columns [birth_pLDDT, death_pLDDT].
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT, N_cc].
        One row per insertion step, in decreasing-pLDDT order.
    """
    plddt = np.asarray(plddt, dtype=np.float64)
    n = len(plddt)
    if n == 0:
        return np.empty((0, 2), dtype=np.float64), np.empty((0, 2), dtype=np.float64)

    # stable sort preserves left-to-right order among ties.
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


# ---------------------------------------------------------------------------
# Size-elder ablation (identical PD to canonical; kept for completeness)
# ---------------------------------------------------------------------------

def build_pd_and_ncc_batch_size_elder(
    plddt: np.ndarray,
    discretize=np.round,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Batch filtration with size-based elder tie-break.

    Identical to ``build_pd_and_ncc`` except: when two components share the
    same birth_pLDDT, the LARGER component (more residues) survives.
    Second tie-break: lower root index.

    In practice produces the same PD as ``build_pd_and_ncc`` because
    within-batch merges always yield birth = death = L regardless of which
    root survives.
    """
    plddt = np.asarray(plddt, dtype=np.float64)
    n = len(plddt)
    if n == 0:
        return np.empty((0, 2), dtype=np.float64), np.empty((0, 2), dtype=np.float64)

    plddt_disc = discretize(plddt).astype(np.int64)
    levels = np.unique(plddt_disc)[::-1]

    parent = np.arange(n, dtype=np.intp)
    birth_of_root = np.zeros(n, dtype=np.float64)
    size = np.ones(n, dtype=np.int64)
    inserted = np.zeros(n, dtype=bool)

    pd_pairs: list[tuple[float, float]] = []
    ncc_rows: list[tuple[float, int]] = []
    n_cc = 0

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)

    def _merge_size(a: int, b: int, death: float) -> None:
        nonlocal n_cc
        ra, rb = find(a), find(b)
        if ra == rb:
            pd_pairs.append((death, death))
            return
        ba, bb = birth_of_root[ra], birth_of_root[rb]
        ra_survives = (
            ba > bb
            or (ba == bb and size[ra] > size[rb])
            or (ba == bb and size[ra] == size[rb] and ra < rb)
        )
        if ra_survives:
            size[ra] += size[rb]
            parent[rb] = ra
            pd_pairs.append((bb, death))
        else:
            size[rb] += size[ra]
            parent[ra] = rb
            pd_pairs.append((ba, death))
        n_cc -= 1

    for level in levels:
        L = float(level)
        batch = np.where(plddt_disc == level)[0]

        for idx in batch:
            i = int(idx)
            parent[i] = i
            birth_of_root[i] = L
            size[i] = 1
            inserted[i] = True
            n_cc += 1

        for idx in batch:
            i = int(idx)
            if i > 0 and inserted[i - 1]:
                _merge_size(i, i - 1, L)
            if i < n - 1 and inserted[i + 1] and plddt_disc[i + 1] > level:
                _merge_size(i, i + 1, L)

        for _ in batch:
            ncc_rows.append((L, n_cc))

    pd_arr = np.array(pd_pairs, dtype=np.float64).reshape(-1, 2)
    return pd_arr, np.array(ncc_rows, dtype=np.float64)
