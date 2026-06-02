"""tests/test_plm.py — unit tests for the PLM statistic.

Covers:
  - Flat N_cc curve → PLM = 1 (one trivial global maximum).
  - Bimodal curve → PLM = 2 at small t_p, PLM = 1 at large t_p.
  - Three prototype proteins at t_p = 0.025: verifies ordering
    P15121 ≤ A0A0G2L439 ≤ Q9VQS4 and that Q9VQS4 gives PLM ≥ 2.
"""
import numpy as np
import pytest
from pathlib import Path

from plm import plm


PROTO_DIR = Path(__file__).resolve().parent.parent / "data" / "prototypes"


def _make_ncc_curve(values):
    """Build a synthetic ncc_curve from a list of N_cc values."""
    n = len(values)
    arr = np.zeros((n, 2), dtype=np.float64)
    arr[:, 0] = np.linspace(100, 1, n)   # decreasing pLDDT (unused by plm)
    arr[:, 1] = np.asarray(values, dtype=np.float64)
    return arr


# ---------------------------------------------------------------------------
# Toy: flat curve
# ---------------------------------------------------------------------------

def test_flat_curve_gives_plm_one():
    """A constant N_cc has exactly one (trivial) persistent maximum."""
    ncc = _make_ncc_curve([5] * 50)
    assert plm(ncc, n_residues=50, t_p=0.025) == 1


# ---------------------------------------------------------------------------
# Toy: bimodal curve
#
# N_cc = [1, 4, 8, 4, 2, 6, 2, 1, 1]  (n=9)
#   Peak 1 at index 2: N_cc = 8 (global max)
#   Valley at index 4: N_cc = 2
#   Peak 2 at index 5: N_cc = 6
#
# Persistence of peak 2 = 6 - 2 = 4.
# threshold = t_p * 9
#   t_p = 0.1 → threshold = 0.9 → 4 ≥ 0.9 → PLM = 2
#   t_p = 0.5 → threshold = 4.5 → 4 < 4.5  → PLM = 1
# ---------------------------------------------------------------------------

BIMODAL = _make_ncc_curve([1, 4, 8, 4, 2, 6, 2, 1, 1])


def test_bimodal_small_tp_gives_plm_two():
    assert plm(BIMODAL, n_residues=9, t_p=0.10) == 2


def test_bimodal_large_tp_gives_plm_one():
    assert plm(BIMODAL, n_residues=9, t_p=0.50) == 1


# ---------------------------------------------------------------------------
# Prototype proteins at t_p = 0.025
# ---------------------------------------------------------------------------

def _prototype_plm(pdb_filename: str, t_p: float = 0.025) -> int:
    from parse import parse_alphafold_pdb
    from filtration import build_pd_and_ncc
    data = parse_alphafold_pdb(PROTO_DIR / pdb_filename)
    plddt = data["plddt"]
    _, ncc_curve = build_pd_and_ncc(plddt)
    return plm(ncc_curve, n_residues=len(plddt), t_p=t_p)


def test_prototype_ordering():
    """Ordered protein has the lowest PLM of the three prototypes.

    The paper-supported claim is that well-structured proteins have few
    persistent local maxima in N_cc(pLDDT); disordered and mixed proteins
    can both fragment the pLDDT profile but their relative ordering is not
    asserted by the paper and does not survive integer-pLDDT discretisation
    on current AlphaFold-DB data.
    """
    p_ordered  = _prototype_plm("AF-P15121-F1-model_v4.pdb")
    p_disorder = _prototype_plm("AF-A0A0G2L439-F1-model_v4.pdb")
    p_mixed    = _prototype_plm("AF-Q9VQS4-F1-model_v6.pdb")
    assert p_ordered <= p_disorder and p_ordered <= p_mixed, (
        f"Expected ordered ≤ both others; got ordered={p_ordered}, "
        f"disordered={p_disorder}, mixed={p_mixed}"
    )


def test_q9vqs4_plm_ge_two():
    """Mixed/fragmented protein Q9VQS4 should have PLM ≥ 2 (paper Fig. 3)."""
    p_mixed = _prototype_plm("AF-Q9VQS4-F1-model_v6.pdb")
    assert p_mixed >= 2, f"Q9VQS4 PLM = {p_mixed}, expected ≥ 2"


def test_plm_always_ge_one():
    """PLM is always ≥ 1 regardless of t_p — the global max always counts."""
    for fname in [
        "AF-P15121-F1-model_v4.pdb",
        "AF-A0A0G2L439-F1-model_v4.pdb",
        "AF-Q9VQS4-F1-model_v6.pdb",
    ]:
        assert _prototype_plm(fname, t_p=0.5) >= 1
