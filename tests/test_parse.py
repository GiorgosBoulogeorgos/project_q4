"""tests/test_parse.py — smoke tests for parse_alphafold_pdb on the three prototypes."""
from pathlib import Path
import numpy as np
import pytest

from parse import parse_alphafold_pdb

PROTOTYPES_DIR = Path(__file__).parent.parent / "data" / "prototypes"

PROTOTYPES = [
    ("P15121",     316, "Ordered"),
    ("A0A0G2L439", 449, "Disordered"),
    ("Q9VQS4",     781, "Mixed"),
]


def _find_pdb(uid: str) -> Path:
    matches = sorted(PROTOTYPES_DIR.glob(f"AF-{uid}-F1-model_v*.pdb"))
    if not matches:
        pytest.skip(f"Prototype PDB for {uid} not found in {PROTOTYPES_DIR}")
    return matches[-1]


@pytest.mark.parametrize("uid,expected_n,label", PROTOTYPES)
def test_parse_returns_correct_length(uid, expected_n, label):
    path = _find_pdb(uid)
    result = parse_alphafold_pdb(path)
    assert result["uniprot_id"] == uid, f"{label}: wrong UniProt ID"
    assert len(result["sequence"]) == expected_n, (
        f"{label}: sequence length {len(result['sequence'])} != {expected_n}"
    )
    assert len(result["plddt"]) == expected_n, (
        f"{label}: pLDDT length {len(result['plddt'])} != {expected_n}"
    )
    assert result["ca_coords"].shape == (expected_n, 3), (
        f"{label}: ca_coords shape {result['ca_coords'].shape}"
    )


@pytest.mark.parametrize("uid,expected_n,label", PROTOTYPES)
def test_plddt_in_range(uid, expected_n, label):
    path = _find_pdb(uid)
    result = parse_alphafold_pdb(path)
    plddt = result["plddt"]
    assert np.all(plddt >= 0.0) and np.all(plddt <= 100.0), (
        f"{label}: pLDDT out of [0, 100] range "
        f"(min={plddt.min():.2f}, max={plddt.max():.2f})"
    )


@pytest.mark.parametrize("uid,expected_n,label", PROTOTYPES)
def test_sequence_is_single_letter_codes(uid, expected_n, label):
    path = _find_pdb(uid)
    result = parse_alphafold_pdb(path)
    valid = set("ACDEFGHIKLMNPQRSTVWYUOX")
    bad = set(result["sequence"]) - valid
    assert not bad, f"{label}: unexpected residue codes {bad}"
