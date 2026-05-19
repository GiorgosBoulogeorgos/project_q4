"""tests/test_statistics.py — unit tests for the four pLDDT statistics."""
import numpy as np
import pytest

from filtration import build_pd_and_ncc
from statistics import f_cp_plus, mean_persistence, persistence_entropy, ncc_max


# ---------------------------------------------------------------------------
# Trivial / analytical checks
# ---------------------------------------------------------------------------

def _make_pd(pairs):
    """Build a PD from a list of (birth, death) tuples."""
    return np.array(pairs, dtype=np.float64)


def test_f_cp_plus_all_positive():
    pd = _make_pd([(80, 60), (70, 50), (90, 40)])
    assert f_cp_plus(pd) == pytest.approx(1.0)


def test_f_cp_plus_all_accretions():
    pd = _make_pd([(50, 50), (30, 30), (20, 20)])
    assert f_cp_plus(pd) == pytest.approx(0.0)


def test_f_cp_plus_mixed():
    pd = _make_pd([(80, 60), (50, 50), (30, 30)])
    assert f_cp_plus(pd) == pytest.approx(1 / 3)


def test_mean_persistence_zero():
    pd = _make_pd([(50, 50), (30, 30)])
    assert mean_persistence(pd) == pytest.approx(0.0)


def test_mean_persistence_known():
    pd = _make_pd([(80, 60), (70, 50)])  # persistences: 20, 20
    assert mean_persistence(pd) == pytest.approx(20.0)


def test_persistence_entropy_uniform():
    """All distinct → maximum entropy → H_p = 1."""
    pd = _make_pd([(80, 60), (70, 50), (90, 30)])  # persistences: 20, 20, 60
    # Two distinct values (20 and 60).
    # P[20] = 2/3, P[60] = 1/3
    # H_raw = -(2/3 * ln(2/3) + 1/3 * ln(1/3))
    # H_p = H_raw / ln(2)
    probs = np.array([2 / 3, 1 / 3])
    expected = float(-np.sum(probs * np.log(probs)) / np.log(2))
    assert persistence_entropy(pd) == pytest.approx(expected, abs=1e-10)


def test_persistence_entropy_all_same():
    """Single distinct value → H_p = 0."""
    pd = _make_pd([(80, 60), (70, 50)])  # both persistence = 20
    assert persistence_entropy(pd) == pytest.approx(0.0)


def test_ncc_max_value():
    ncc_curve = np.array([(90, 1), (80, 2), (70, 3), (60, 2), (50, 1)], dtype=np.float64)
    assert ncc_max(ncc_curve) == 3


# ---------------------------------------------------------------------------
# Null model: H_p ≈ 0.45  (paper, n = 1000)
# ---------------------------------------------------------------------------

def test_null_model_persistence_entropy():
    """
    For n=1000 i.i.d. Uniform[0,100] pLDDT, H_p should converge to ≈ 0.45.
    We use a single large trial (n=2000) to reduce variance.
    """
    rng = np.random.default_rng(7)
    plddt = rng.uniform(0, 100, size=2000)
    pd, _ = build_pd_and_ncc(plddt)
    hp = persistence_entropy(pd)
    assert abs(hp - 0.45) < 0.05, (
        f"Null model H_p = {hp:.4f}; expected ≈ 0.45 ± 0.05"
    )


# ---------------------------------------------------------------------------
# Round-trip with filtration on a tiny hand-crafted example
# ---------------------------------------------------------------------------

def test_two_node_path():
    """Two residues: one PD pair (accretion), f+_cp=0, p_bar=0, H_p=0, ncc_max=1."""
    plddt = np.array([80.0, 60.0])
    pd, ncc = build_pd_and_ncc(plddt)
    assert f_cp_plus(pd) == pytest.approx(0.0)
    assert mean_persistence(pd) == pytest.approx(0.0)
    assert persistence_entropy(pd) == pytest.approx(0.0)
    assert ncc_max(ncc) == 1
