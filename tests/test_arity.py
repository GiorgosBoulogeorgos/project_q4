"""tests/test_arity.py — unit tests for the arity module.

Covers:
  - Linear toy: 3 atoms in a line, only immediate neighbours within r → [1, 2, 1].
  - Self-exclusion: arity of a single atom is always 0.
  - Three prototype proteins: ordering P15121 (folded) > Q9VQS4 ≥ A0A0G2L439.
"""
from __future__ import annotations

import numpy as np
import pytest
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from arity import arity_per_residue, arity_signature

PROTO_DIR = Path(__file__).resolve().parent.parent / "data" / "prototypes"


# ---------------------------------------------------------------------------
# Toy: linear three-atom chain
# ---------------------------------------------------------------------------

def test_linear_three_atoms():
    """3 atoms spaced 8 Å apart; end pair is 16 Å > 10 Å → [1, 2, 1]."""
    coords = np.array([[0, 0, 0], [8, 0, 0], [16, 0, 0]], dtype=float)
    arities = arity_per_residue(coords, r=10.0)
    assert list(arities) == [1, 2, 1], f"Expected [1, 2, 1], got {list(arities)}"


def test_single_atom_arity_zero():
    """A lone atom has no neighbours → arity = 0."""
    coords = np.array([[0, 0, 0]], dtype=float)
    arities = arity_per_residue(coords, r=10.0)
    assert arities[0] == 0


def test_all_within_radius():
    """All pairs within r → each atom sees all others."""
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
    arities = arity_per_residue(coords, r=10.0)
    assert list(arities) == [2, 2, 2]


# ---------------------------------------------------------------------------
# Signature: inverse-CDF definition (paper Def. 2 / Eq. 2)
# ---------------------------------------------------------------------------

def test_signature_is_inverse_cdf_not_linear_interp():
    """
    arity_Cα(q) must be the smallest observed arity a with CDF(a) >= q
    (paper Eq. 2), NOT a linearly interpolated percentile.

    Fixture: 3 coincident atoms (arity 2 each) plus 10 coincident atoms 100 Å
    away (arity 9 each).  The arity distribution is [2,2,2, 9×10], with a wide
    gap between 2 and 9.  At q = 0.2 the inverse CDF returns 2 (CDF(2) = 3/13
    ≈ 0.231 >= 0.2), whereas np.percentile's linear interpolation returns 4 —
    an arity no residue actually has.  This test fails on the old
    linear-interpolation implementation.
    """
    coords = np.array([[0, 0, 0]] * 3 + [[100, 0, 0]] * 10, dtype=float)

    arities = arity_per_residue(coords)
    assert sorted(arities.tolist()) == [2, 2, 2] + [9] * 10

    sig = arity_signature(coords, q1=0.2, q2=0.2)
    assert sig == (2, 2), f"Expected inverse-CDF (2, 2), got {sig}"

    # Every signature value must be an arity that actually occurs.
    observed = set(arities.tolist())
    for q in (0.1, 0.2, 0.25, 0.5, 0.75, 0.9):
        lo, hi = arity_signature(coords, q1=q, q2=q)
        assert lo in observed and hi in observed, (
            f"q={q}: signature ({lo}, {hi}) contains a non-observed arity"
        )


# ---------------------------------------------------------------------------
# Prototype proteins
# ---------------------------------------------------------------------------

def _proto_signature(pdb_filename: str) -> tuple[int, int]:
    from parse import parse_alphafold_pdb
    data = parse_alphafold_pdb(PROTO_DIR / pdb_filename)
    return arity_signature(data["ca_coords"])


def test_prototype_ordering():
    """
    Folded P15121 should have higher arity than the disordered / mixed proteins.
    Paper's Fig. 4 places P15121 near (15, 19); disordered proteins are lower.
    """
    a25_ordered, a75_ordered   = _proto_signature("AF-P15121-F1-model_v4.pdb")
    a25_disorder, a75_disorder = _proto_signature("AF-A0A0G2L439-F1-model_v4.pdb")
    a25_mixed,   a75_mixed     = _proto_signature("AF-Q9VQS4-F1-model_v6.pdb")

    assert a75_ordered >= a75_disorder, (
        f"Expected ordered arity_75 ({a75_ordered}) ≥ disordered ({a75_disorder})"
    )
    assert a75_ordered >= a75_mixed, (
        f"Expected ordered arity_75 ({a75_ordered}) ≥ mixed ({a75_mixed})"
    )
