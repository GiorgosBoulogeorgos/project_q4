"""tests/test_filtration_batch.py — ablation tests for batch variants.

Tests build_pd_and_ncc with discretize=np.floor (Experiment C) and
build_pd_and_ncc_batch_size_elder (Experiment B ablation).

Core invariant tests for the canonical batch function live in test_filtration.py.
"""
import numpy as np
import pytest

from filtration import build_pd_and_ncc, build_pd_and_ncc_batch_size_elder


def test_floor_pd_shape():
    rng = np.random.default_rng(3)
    plddt = rng.uniform(0, 100, size=50)
    pd, ncc = build_pd_and_ncc(plddt, discretize=np.floor)
    assert pd.shape == (49, 2)
    assert ncc.shape == (50, 2)
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_floor_birth_ge_death():
    rng = np.random.default_rng(8)
    pd, _ = build_pd_and_ncc(rng.uniform(0, 100, size=300), discretize=np.floor)
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_size_elder_pd_shape():
    rng = np.random.default_rng(9)
    plddt = rng.uniform(0, 100, size=50)
    pd, ncc = build_pd_and_ncc_batch_size_elder(plddt)
    assert pd.shape == (49, 2)
    assert ncc.shape == (50, 2)


def test_size_elder_birth_ge_death():
    rng = np.random.default_rng(10)
    pd, _ = build_pd_and_ncc_batch_size_elder(rng.uniform(0, 100, size=300))
    assert np.all(pd[:, 0] >= pd[:, 1])


def test_size_elder_determinism():
    rng = np.random.default_rng(7)
    plddt = rng.uniform(0, 100, size=200)
    pd1, _ = build_pd_and_ncc_batch_size_elder(plddt)
    pd2, _ = build_pd_and_ncc_batch_size_elder(plddt)
    assert np.array_equal(pd1, pd2)
