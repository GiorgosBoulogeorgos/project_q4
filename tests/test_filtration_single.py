"""tests/test_filtration_single.py — unit tests for build_pd_and_ncc_single.

This is the original single-insertion algorithm (one residue per step, raw
float pLDDT).  Kept as a reference/ablation test suite.  The canonical
algorithm is tested in test_filtration.py.
"""
import numpy as np
import pytest

from filtration import build_pd_and_ncc_single


TOY_PLDDT = np.array([10, 50, 30, 90, 20, 80, 40, 70, 60, 100], dtype=np.float64)


def test_toy_pd_size():
    pd, ncc = build_pd_and_ncc_single(TOY_PLDDT)
    assert pd.shape == (len(TOY_PLDDT) - 1, 2)


def test_toy_ncc_size():
    pd, ncc = build_pd_and_ncc_single(TOY_PLDDT)
    assert ncc.shape == (len(TOY_PLDDT), 2)


def test_toy_birth_ge_death():
    pd, _ = build_pd_and_ncc_single(TOY_PLDDT)
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_toy_ncc_positive():
    _, ncc = build_pd_and_ncc_single(TOY_PLDDT)
    assert np.all(ncc[:, 1] >= 1)


def test_toy_ncc_starts_at_one():
    _, ncc = build_pd_and_ncc_single(TOY_PLDDT)
    assert int(ncc[0, 1]) == 1


def test_toy_final_ncc_is_one():
    _, ncc = build_pd_and_ncc_single(TOY_PLDDT)
    assert int(ncc[-1, 1]) == 1


def test_positive_persistence_count():
    pd, _ = build_pd_and_ncc_single(TOY_PLDDT)
    assert int(np.sum(pd[:, 0] > pd[:, 1])) >= 1


def test_determinism():
    rng = np.random.default_rng(0)
    plddt = rng.uniform(0, 100, size=200)
    pd1, ncc1 = build_pd_and_ncc_single(plddt)
    pd2, ncc2 = build_pd_and_ncc_single(plddt)
    assert np.array_equal(pd1, pd2)
    assert np.array_equal(ncc1, ncc2)


def test_null_model_f_cp_plus():
    from statistics import f_cp_plus

    rng = np.random.default_rng(42)
    results = []
    for _ in range(20):
        plddt = rng.uniform(0, 100, size=1000)
        pd, _ = build_pd_and_ncc_single(plddt)
        results.append(f_cp_plus(pd))

    mean_f = float(np.mean(results))
    assert abs(mean_f - 1 / 3) < 0.05


def test_single_residue():
    pd, ncc = build_pd_and_ncc_single(np.array([75.0]))
    assert pd.shape == (0, 2)
    assert ncc.shape == (1, 2)
    assert int(ncc[0, 1]) == 1


def test_two_residues():
    pd, ncc = build_pd_and_ncc_single(np.array([80.0, 60.0]))
    assert pd.shape == (1, 2)
    assert pd[0, 0] >= pd[0, 1]
    assert int(ncc[-1, 1]) == 1
