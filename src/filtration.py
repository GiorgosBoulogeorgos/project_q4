"""
filtration.py — path-graph pLDDT filtration with Union-Find (Elder Rule).

A polypeptide of n residues is modelled as a path graph
0 — 1 — 2 — … — (n−1); each vertex carries a pLDDT score.  The filtration
inserts residues in DECREASING pLDDT order; edges are activated when both
endpoints are present.  This is the standard sublevel-set filtration described
in Cazals & Sarti (2025), using u = −pLDDT as the filter function.

pLDDT discretisation
--------------------
By default, ``build_pd_and_ncc`` rounds pLDDT to integers in [0, 100] before
running the filtration.  This is what AlphaFold itself reports internally and
is the form of pLDDT for which the paper's null-model H_p values match across
all three sample sizes (n=100, 1000, 10 000); the raw-float baseline only
matches at n=1000 by coincidence (see ``docs/hp_investigation.md`` and the
``project_hp_discretisation`` memory).  The flag ``discretise=False`` is
provided so the original raw-float baseline remains accessible for
side-by-side comparisons.

Batched insertion
-----------------
By default, residues sharing the same integer pLDDT level are inserted
\\emph{in batch}: every singleton is created first, then every incident
edge is activated exactly once.  This matches the paper's description that
"pLDDT values come in batches, so that local insertions into stretches
along the sequence prevent the maximum number of local maxima to rise"
(Cazals & Sarti 2025, §3.1).  The resulting persistence diagram is
mathematically identical to the sequential ascending-index variant under
the Elder Rule, but the N_cc curve has flat plateaus instead of intra-batch
transients — which significantly reduces PLM counts.  Set
``batched=False`` to recover the sequential one-residue-per-step
behaviour for ablation.

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
    One row per residue: (pLDDT_value, N_cc).  Under ``batched=True``
    every residue in a batch carries the post-batch N_cc value, producing
    a piecewise-constant curve; under ``batched=False`` the curve is
    recorded after each individual insertion, exposing intra-batch
    transients.  Rows are ordered by decreasing pLDDT.
"""

from __future__ import annotations

import numpy as np


def build_pd_and_ncc(
    plddt: np.ndarray,
    discretise: bool = True,
    batched: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the path-graph pLDDT filtration.

    Parameters
    ----------
    plddt:
        1-D array of per-residue pLDDT values (0–100 scale), in residue order.
    discretise:
        When True (default), round pLDDT to integers on [0, 100] before the
        filtration.  This matches the discretisation that AlphaFold itself
        reports and is required to reproduce the paper's null-model H_p
        values at n=10 000.  Set False to use raw float pLDDT for ablation
        or sensitivity studies.
    batched:
        When True (default), residues sharing the same pLDDT level are
        inserted in batch: every singleton is created first, then every
        incident edge is activated exactly once in canonical (low, high)
        order.  Matches the paper's "pLDDT values come in batches"
        statement.  When False, residues are processed one at a time in
        decreasing-pLDDT order with ties broken by ascending index
        (stable sort) — the original sequential variant.

    Returns
    -------
    pd:
        Array of shape (n−1, 2) with columns [birth_pLDDT, death_pLDDT].
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT_value, N_cc], rows
        ordered by decreasing pLDDT.
    """
    plddt = np.asarray(plddt, dtype=np.float64)
    if discretise:
        plddt = np.round(plddt)
    n = len(plddt)
    if n == 0:
        return np.empty((0, 2), dtype=np.float64), np.empty((0, 2), dtype=np.float64)

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

    if batched:
        # ── batched insertion: process by descending pLDDT level ─────────
        for v in sorted(set(plddt.tolist()), reverse=True):
            idxs = sorted(int(j) for j in np.where(plddt == v)[0])
            # (1) introduce every singleton in this level
            for i in idxs:
                parent[i] = i
                birth_of_root[i] = v
                inserted[i] = True
                n_cc += 1
            # (2) collect every incident inserted edge as a canonical
            #     (low, high) pair so each edge fires at most once
            edges: set[tuple[int, int]] = set()
            for i in idxs:
                if i - 1 >= 0 and inserted[i - 1]:
                    edges.add((i - 1, i))
                if i + 1 < n and inserted[i + 1]:
                    edges.add((i, i + 1))
            # (3) merge along each edge once, in canonical order
            for a, b in sorted(edges):
                _merge(a, b, v)
            # (4) record n_cc for every residue in this batch
            for i in idxs:
                ncc_curve.append((v, n_cc))
    else:
        # ── sequential insertion, ascending index within a tie ──────────
        order = np.argsort(-plddt, kind="stable")
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
