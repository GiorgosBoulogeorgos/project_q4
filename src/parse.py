"""
parse.py — AlphaFold PDB/mmCIF → per-residue pLDDT array and Cα coordinates.

pLDDT is stored in the B-factor column (0–100 scale) of AlphaFold-DB files.
Cα coordinates are extracted for the Q1 arity-map extension.
"""

from __future__ import annotations

import re
from pathlib import Path

import gemmi
import numpy as np

# ---------------------------------------------------------------------------
# Amino-acid look-up table (standard 20 + selenocysteine/pyrrolysine)
# ---------------------------------------------------------------------------

_THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",  # rare but valid
}

_AF_NAME_RE = re.compile(
    r"AF-(?P<uniprot>[^-]+)-F\d+-model_v\d+\.(pdb|cif)$", re.IGNORECASE
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_alphafold_pdb(path: Path) -> dict:
    """
    Parse an AlphaFold-DB PDB (or mmCIF) file.

    Parameters
    ----------
    path:
        Path to the ``.pdb`` or ``.cif`` file.

    Returns
    -------
    dict with keys:
        ``uniprot_id`` (str)
            Accession parsed from the filename, e.g. ``"P15121"``.
        ``sequence`` (str)
            One-letter amino-acid sequence; non-standard residues become
            ``"X"``.
        ``plddt`` (np.ndarray, shape (n,), float64)
            Per-residue pLDDT scores in the [0, 100] range, taken from
            the B-factor column.  One value per Cα residue.
        ``ca_coords`` (np.ndarray, shape (n, 3), float64)
            Cα (x, y, z) coordinates in Å, in chain/residue order.
    """
    path = Path(path)

    m = _AF_NAME_RE.match(path.name)
    uniprot_id: str = m.group("uniprot") if m else path.stem

    st = gemmi.read_structure(str(path))
    model = st[0]

    seq_chars: list[str] = []
    plddt_vals: list[float] = []
    coords: list[list[float]] = []

    for chain in model:
        for res in chain:
            ca = res.find_atom("CA", "\0")
            if ca is None:
                continue
            seq_chars.append(_THREE_TO_ONE.get(res.name, "X"))
            plddt_vals.append(float(ca.b_iso))
            coords.append([ca.pos.x, ca.pos.y, ca.pos.z])

    return {
        "uniprot_id": uniprot_id,
        "sequence": "".join(seq_chars),
        "plddt": np.array(plddt_vals, dtype=np.float64),
        "ca_coords": np.array(coords, dtype=np.float64),
    }
