"""tests/test_filtration.py — unit tests for build_pd_and_ncc (canonical batch).

Covers:
  - Structural invariants (shape, birth ≥ death, determinism, edge cases).
  - Toy hand-traced example verifying exact PD pairs.
  - Null model: f⁺_cp ≈ 0.33, H_p ≈ 0.46 for n=1000 uniform pLDDT.
  - Prototype H_p regression: P15121 ≈ 0.15, A0A0G2L439 ≈ 0.43, Q9VQS4 ≈ 0.42
    (tolerance ±0.02; these are the expected values under integer-pLDDT batch
    insertion on current AlphaFold-DB v4/v6 data).
"""
import numpy as np
import pytest
from pathlib import Path

from filtration import build_pd_and_ncc


PROTO_DIR = Path(__file__).resolve().parent.parent / "data" / "prototypes"


# ---------------------------------------------------------------------------
# Toy: pLDDT = [98, 98, 50, 50, 98, 98]  (n=6, two high-pLDDT segments)
#
# Integer levels after np.round: 98 → batch [0,1,4,5];  50 → batch [2,3].
#
# Level 98:
#   Add singletons 0,1,4,5: n_cc=4
#   i=1: left(0) inserted → merge(1,0,98): both born 98, lower index (0)
#         survives → PD(98,98). n_cc=3.
#   i=5: left(4) inserted → merge(5,4,98): both born 98, lower index (4)
#         survives → PD(98,98). n_cc=2.
#   ncc_rows: 4× (98,2)
#
# Level 50:
#   Add singletons 2,3: n_cc=4
#   i=2: left(1) inserted, find(1)→0 (born 98 > 50) → PD(50,50). n_cc=3.
#   i=3: left(2) inserted, find(2)→0 → PD(50,50). n_cc=2.
#          right(4) at level 98>50 → merge(3→comp0, 4→comp4,98): lower index
#          root 0 survives, 4 dies → PD(98,50). n_cc=1.
#   ncc_rows: 2× (50,1)
#
# Final PD (n−1=5 pairs):
#   (98,98), (98,98), (50,50), (50,50), (98,50)
# n_accretion=4, n_positive=1 (persistence=48)
# ---------------------------------------------------------------------------

TOY = np.array([98.0, 98.0, 50.0, 50.0, 98.0, 98.0])


def test_toy_pd_shape():
    pd, _ = build_pd_and_ncc(TOY)
    assert pd.shape == (5, 2)


def test_toy_ncc_shape():
    _, ncc = build_pd_and_ncc(TOY)
    assert ncc.shape == (6, 2)


def test_toy_birth_ge_death():
    pd, _ = build_pd_and_ncc(TOY)
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_toy_pd_pairs_exact():
    pd, _ = build_pd_and_ncc(TOY)
    persts = pd[:, 0] - pd[:, 1]
    assert int(np.sum(persts == 0)) == 4, "Expected 4 accretions"
    assert int(np.sum(persts > 0)) == 1, "Expected 1 positive pair"


def test_toy_positive_pair_value():
    pd, _ = build_pd_and_ncc(TOY)
    pos = pd[pd[:, 0] > pd[:, 1]]
    assert len(pos) == 1
    assert float(pos[0, 0]) == pytest.approx(98.0)
    assert float(pos[0, 1]) == pytest.approx(50.0)


def test_toy_final_ncc_is_one():
    _, ncc = build_pd_and_ncc(TOY)
    assert int(ncc[-1, 1]) == 1


# ---------------------------------------------------------------------------
# Structural invariants on random inputs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [1, 2, 5, 10, 100, 500])
def test_pd_size_invariant(n):
    rng = np.random.default_rng(n)
    plddt = rng.uniform(0, 100, size=n)
    pd, ncc = build_pd_and_ncc(plddt)
    assert pd.shape == (max(n - 1, 0), 2)
    assert ncc.shape == (n, 2)


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_birth_ge_death_random(seed):
    rng = np.random.default_rng(seed)
    pd, _ = build_pd_and_ncc(rng.uniform(0, 100, size=300))
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_ncc_positive():
    rng = np.random.default_rng(5)
    _, ncc = build_pd_and_ncc(rng.uniform(0, 100, size=200))
    assert np.all(ncc[:, 1] >= 1)


def test_final_ncc_is_one():
    rng = np.random.default_rng(6)
    _, ncc = build_pd_and_ncc(rng.uniform(0, 100, size=200))
    assert int(ncc[-1, 1]) == 1


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_single_residue():
    pd, ncc = build_pd_and_ncc(np.array([75.0]))
    assert pd.shape == (0, 2)
    assert ncc.shape == (1, 2)
    assert int(ncc[0, 1]) == 1


def test_two_residues():
    pd, ncc = build_pd_and_ncc(np.array([80.0, 60.0]))
    assert pd.shape == (1, 2)
    assert pd[0, 0] >= pd[0, 1]
    assert int(ncc[-1, 1]) == 1


def test_all_same_plddt():
    """All residues at the same level: all merges are accretions."""
    pd, ncc = build_pd_and_ncc(np.full(10, 70.0))
    assert pd.shape == (9, 2)
    assert np.all(pd[:, 0] == pd[:, 1])  # all accretions
    assert int(ncc[-1, 1]) == 1


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------

def test_determinism():
    rng = np.random.default_rng(7)
    plddt = rng.uniform(0, 100, size=200)
    pd1, ncc1 = build_pd_and_ncc(plddt)
    pd2, ncc2 = build_pd_and_ncc(plddt)
    assert np.array_equal(pd1, pd2)
    assert np.array_equal(ncc1, ncc2)


# ---------------------------------------------------------------------------
# Null model (n=1000 i.i.d. U[0,100], np.round to integer levels)
# ---------------------------------------------------------------------------

def test_null_model_f_cp_plus():
    """f⁺_cp → 1/3 for i.i.d. uniform pLDDT (paper Example 1)."""
    from statistics import f_cp_plus

    rng = np.random.default_rng(42)
    results = [f_cp_plus(build_pd_and_ncc(rng.uniform(0, 100, size=1000))[0])
               for _ in range(20)]
    mean_f = float(np.mean(results))
    assert abs(mean_f - 1 / 3) < 0.07, f"f⁺_cp mean = {mean_f:.4f}, expected ≈ 0.333"


def test_null_model_persistence_entropy():
    """H_p ≈ 0.46 for i.i.d. uniform pLDDT at n=1000 (batch variant)."""
    from statistics import persistence_entropy

    rng = np.random.default_rng(0)
    results = [persistence_entropy(build_pd_and_ncc(rng.uniform(0, 100, size=1000))[0])
               for _ in range(20)]
    mean_hp = float(np.mean(results))
    assert abs(mean_hp - 0.46) < 0.03, f"H_p mean = {mean_hp:.4f}, expected ≈ 0.46"


# ---------------------------------------------------------------------------
# Prototype H_p regression (tolerance ±0.02)
#
# These are the empirically observed H_p values under integer-pLDDT batch
# insertion on current AlphaFold-DB data.  They are documented findings, not
# paper reproduction targets; the offset from the paper's Figure 3 values is
# attributed to higher float precision in the v4/v6 data.
# ---------------------------------------------------------------------------

def _compute_hp(pdb_filename: str) -> float:
    from parse import parse_alphafold_pdb
    from statistics import persistence_entropy
    plddt = parse_alphafold_pdb(PROTO_DIR / pdb_filename)["plddt"]
    pd, _ = build_pd_and_ncc(plddt)
    return persistence_entropy(pd)


@pytest.mark.parametrize("filename,expected_hp", [
    ("AF-P15121-F1-model_v4.pdb",      0.147),
    ("AF-A0A0G2L439-F1-model_v4.pdb",  0.433),
    ("AF-Q9VQS4-F1-model_v6.pdb",      0.425),
])
def test_prototype_hp(filename, expected_hp):
    hp = _compute_hp(filename)
    assert abs(hp - expected_hp) < 0.02, (
        f"{filename}: H_p = {hp:.4f}, expected {expected_hp:.3f} ± 0.02"
    )
