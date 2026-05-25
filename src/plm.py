"""
plm.py — Persistent Local Maxima (PLM) of the N_cc curve.

PLM(t_p) counts the number of local maxima of the connected-component count
function N_cc(pLDDT) that survive Morse-Smale simplification at relative
persistence threshold t_p.

Implementation
--------------
1. Extract the N_cc column from the ncc_curve returned by build_pd_and_ncc.
2. Feed the *negated* N_cc array to gudhi's CubicalComplex.  Negation turns
   local maxima of N_cc into local minima of -N_cc, which are the features
   captured by the sublevel-set 0-dimensional persistence diagram.
3. The persistence of each 0-dim pair equals N_cc_peak - N_cc_valley (in
   integer N_cc units).
4. PLM(t_p) = number of pairs with persistence ≥ t_p * n_residues.
   The global maximum (birth = -max(N_cc), death = +∞) always contributes
   one count.

gudhi is the only external TDA library used in this project; all other
persistence computation is implemented from scratch in filtration.py.

Reference: Cazals & Sarti (2025), §4.3 — PLM(t_ν) with t_ν = n · t_p.
"""

from __future__ import annotations

import numpy as np
import gudhi


def plm(
    ncc_curve: np.ndarray,
    n_residues: int,
    t_p: float = 0.025,
) -> int:
    """
    Count persistent local maxima of N_cc after simplification at t_p.

    Parameters
    ----------
    ncc_curve:
        Array of shape (n, 2) with columns [pLDDT, N_cc], as returned by
        build_pd_and_ncc.
    n_residues:
        Number of residues in the protein.  Used to scale the persistence
        threshold: t_ν = t_p * n_residues.
    t_p:
        Relative persistence threshold.  Default 0.025 (paper convention).

    Returns
    -------
    int
        Number of local maxima of N_cc with persistence ≥ t_p * n_residues.
        Always ≥ 1 (the global maximum is always counted).
    """
    ncc_values = ncc_curve[:, 1].astype(np.float64)

    # Feed -N_cc so that local maxima of N_cc become local minima of -N_cc,
    # which are born first in the sublevel-set filtration.
    neg_ncc = -ncc_values

    cc = gudhi.CubicalComplex(top_dimensional_cells=neg_ncc)
    pairs = cc.persistence()  # list of (dim, (birth, death))

    threshold = t_p * n_residues
    count = 0
    for dim, (birth, death) in pairs:
        if dim != 0:
            continue
        # persistence in -N_cc units = death - birth = -(valley) - (-(peak))
        #                             = N_cc_peak - N_cc_valley
        persistence = death - birth  # finite or inf
        if persistence >= threshold:
            count += 1

    return count


def plm_multi(
    ncc_curve: np.ndarray,
    n_residues: int,
    t_ps: tuple[float, ...] = (0.020, 0.025, 0.030),
) -> list[int]:
    """
    Compute PLM for multiple t_p values in a single gudhi call.

    Returns a list of counts in the same order as t_ps.
    """
    ncc_values = ncc_curve[:, 1].astype(np.float64)
    cc = gudhi.CubicalComplex(top_dimensional_cells=-ncc_values)
    pairs = cc.persistence()
    persistences = [death - birth for dim, (birth, death) in pairs if dim == 0]
    return [
        sum(1 for p in persistences if p >= t_p * n_residues)
        for t_p in t_ps
    ]


def plm_peak_pLDDTs(
    ncc_curve: np.ndarray,
    n_residues: int,
    t_p: float = 0.025,
) -> list[float]:
    """
    Return the pLDDT levels of persistent local maxima of N_cc at threshold t_p.

    Useful for annotating the N_cc curve plot.
    """
    ncc_values = ncc_curve[:, 1].astype(np.float64)
    plddt_levels = ncc_curve[:, 0]

    cc = gudhi.CubicalComplex(top_dimensional_cells=-ncc_values)
    pairs = cc.persistence()
    threshold = t_p * n_residues

    peak_pLDDTs = []
    for dim, (birth, death) in pairs:
        if dim != 0:
            continue
        if death - birth < threshold:
            continue
        # birth = -N_cc[peak], so N_cc[peak] = -birth
        peak_ncc = -birth
        # find the index where N_cc equals this value AND is a local maximum
        candidates = np.where(np.isclose(ncc_values, peak_ncc, atol=0.5))[0]
        if len(candidates) == 0:
            continue
        # pick the candidate where N_cc is locally maximal
        best = candidates[np.argmax(ncc_values[candidates])]
        peak_pLDDTs.append(float(plddt_levels[best]))

    return peak_pLDDTs
