"""
statistics.py — per-protein pLDDT fragmentation statistics.

Computes from the persistence diagram (PD) and N_cc curve:
    f+_cp     fraction of PD points with strictly positive persistence
    p_bar     mean persistence (mean of birth − death over all n−1 pairs)
    H_p       normalised Shannon entropy of the persistence distribution (Eq. 1)
    ncc_max   peak connected-component count over the filtration sweep

PLM is implemented separately in plm.py (requires gudhi's CubicalComplex).

Formula reference
-----------------
Let P = set of DISTINCT persistence values p_i = birth_i − death_i.
For each distinct value v ∈ P: P[v] = (# PD pairs with persistence v) / m,
where m = n − 1 is the total number of PD pairs.

  H_p = − Σ_{v ∈ P} P[v] · ln P[v]  /  ln |P|       (Eq. 1, Cazals & Sarti)

Edge case: if |P| == 1 (all pairs have identical persistence), H_p = 0.
"""

from __future__ import annotations

import numpy as np


def f_cp_plus(pd: np.ndarray) -> float:
    """
    Fraction of PD pairs with strictly positive persistence.

    f+_cp = m' / m, where m' = #{pairs with birth > death}, m = n−1.
    Accretions (birth == death) contribute to m but not to m'.

    Parameters
    ----------
    pd:
        Persistence diagram, shape (m, 2): columns [birth, death].

    Returns
    -------
    float in [0, 1].
    """
    if len(pd) == 0:
        return 0.0
    return float(np.sum(pd[:, 0] > pd[:, 1])) / len(pd)


def mean_persistence(pd: np.ndarray) -> float:
    """
    Mean persistence p̄ = (1/m) Σ (birth_i − death_i).

    Parameters
    ----------
    pd:
        Persistence diagram, shape (m, 2).

    Returns
    -------
    float >= 0.
    """
    if len(pd) == 0:
        return 0.0
    return float(np.mean(pd[:, 0] - pd[:, 1]))


def persistence_entropy(pd: np.ndarray) -> float:
    """
    Normalised Shannon entropy H_p of the persistence distribution (Eq. 1).

    H_p = − Σ_{v ∈ P} P[v] · ln P[v]  /  ln |P|

    where P is the set of DISTINCT persistence values and P[v] = count(v)/m.

    Parameters
    ----------
    pd:
        Persistence diagram, shape (m, 2).

    Returns
    -------
    float in [0, 1].  Returns 0.0 when m == 0 or |P| == 1.
    """
    if len(pd) == 0:
        return 0.0

    persistences = pd[:, 0] - pd[:, 1]   # shape (m,)
    m = len(persistences)

    _, counts = np.unique(persistences, return_counts=True)
    n_distinct = len(counts)

    if n_distinct <= 1:
        return 0.0

    probs = counts / m
    raw_entropy = -float(np.sum(probs * np.log(probs)))
    return raw_entropy / np.log(n_distinct)


def ncc_max(ncc_curve: np.ndarray) -> int:
    """
    Peak connected-component count N_cc^max over the filtration sweep.

    Parameters
    ----------
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT, N_cc].

    Returns
    -------
    int.
    """
    if len(ncc_curve) == 0:
        return 0
    return int(np.max(ncc_curve[:, 1]))
