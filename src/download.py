"""
download.py — fetch AlphaFold-DB structures and proteome tarballs.

AlphaFold-DB v4 URL patterns
-----------------------------
Single structure (PDB):
  https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F{FRAG}-model_v4.pdb

Full proteome tarball (EBI FTP):
  https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
      UP000005640_9606_HUMAN_v4.tar   (H. sapiens, ~11 GB)

Usage
-----
  python src/download.py          # smoke-test: downloads the 3 prototype PDBs
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import requests
from tqdm import tqdm

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# REST API — returns JSON with pdbUrl, latestVersion, etc.
_AF_API = "https://alphafold.ebi.ac.uk/api/prediction"

# Static file host (individual PDB files; version discovered via API).
_AF_BASE = "https://alphafold.ebi.ac.uk/files"

# Versions to try when the API has no record of a protein (e.g. non-reference
# proteome entries like zebrafish A0A0G2L439).  Tried in order; first 200 wins.
_FALLBACK_VERSIONS = (6, 5, 4, 3)

# EBI FTP — proteome tarballs live under /v4/, not /latest/.
# The v4 proteome tarball (AlphaFold-DB database v4, August 2025) is used
# to match the paper's scope.  Individual structure files are served at their
# current version (v6 as of 2026-05) via the API.
_PROTEOME_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/v4"

_CHUNK = 1 << 20  # 1 MiB streaming chunk

# Map species code → UniProt proteome identifier (without version suffix).
# Source: https://alphafold.ebi.ac.uk/download  /  EBI FTP /v4/
_PROTEOME_IDS: dict[str, str] = {
    "HUMAN": "UP000005640_9606_HUMAN",
    "MOUSE": "UP000000589_10090_MOUSE",
    "RAT":   "UP000002494_10116_RAT",
    "YEAST": "UP000002311_559292_YEAST",
}

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _resume_get(url: str, dest: Path) -> Path:
    """
    Stream *url* to *dest*, resuming from the current file size if it already
    exists.  Uses the HTTP ``Range`` header so interrupted downloads do not
    restart from zero.

    HTTP status meanings handled here:
      206 Partial Content      → resume succeeded, append remaining bytes.
      200 OK                   → server ignored Range (or file was empty).
      416 Range Not Satisfiable → file is already complete, nothing to do.
    """
    existing = dest.stat().st_size if dest.exists() else 0
    headers: dict[str, str] = {}
    if existing:
        headers["Range"] = f"bytes={existing}-"

    with requests.get(url, headers=headers, stream=True, timeout=60) as resp:
        if resp.status_code == 416:
            log.debug("Already complete: %s", dest.name)
            return dest
        resp.raise_for_status()

        content_length = int(resp.headers.get("Content-Length", 0))
        total = existing + content_length
        mode = "ab" if existing else "wb"

        with (
            open(dest, mode) as fh,
            tqdm(
                total=total or None,
                initial=existing,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=dest.name,
                leave=False,
            ) as bar,
        ):
            for chunk in resp.iter_content(_CHUNK):
                fh.write(chunk)
                bar.update(len(chunk))

    return dest


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def download_structure(
    uniprot_id: str,
    out_dir: Path,
    *,
    fragment: int = 1,
    version: int | None = None,
) -> Path:
    """
    Fetch a single AlphaFold-DB prediction by UniProt accession.

    Parameters
    ----------
    uniprot_id:
        UniProt accession, e.g. ``"P15121"``.
    out_dir:
        Directory where the ``.pdb`` file will be written (created if absent).
    fragment:
        Fragment index (1-based).  AlphaFold-DB v4 splits chains longer than
        2 700 residues into overlapping 1 400-residue fragments.  For most
        proteins ``fragment=1`` is the full chain.
    version:
        Model version to pin (e.g. ``4``).  If ``None`` (default), the
        AlphaFold-DB REST API is queried to discover the current latest
        version.  Note: as of 2026, v4 individual files are no longer hosted
        at the EBI; the API returns v6 URLs.

    Returns
    -------
    Path
        Absolute path to the downloaded ``.pdb`` file.

    Notes
    -----
    URL pattern (discovered via API)::

        https://alphafold.ebi.ac.uk/files/AF-P15121-F1-model_v6.pdb

    API endpoint::

        https://alphafold.ebi.ac.uk/api/prediction/P15121
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Short-circuit: if any version of this file is already cached locally,
    # return it without hitting the network.  This handles proteins extracted
    # from a proteome tarball at a version different from what the API returns
    # (e.g. A0A0G2L439 from the Zebrafish v4 tarball).
    cached = sorted(
        out_dir.glob(f"AF-{uniprot_id.upper()}-F{fragment}-model_v*.pdb")
    )
    if cached:
        log.info("Using cached: %s", cached[-1].name)
        return cached[-1]

    if version is None:
        # Resolve via the REST API; fall back to probing file versions directly
        # for entries not indexed by the API (e.g. non-reference-proteome
        # organisms such as zebrafish A0A0G2L439).
        api_resp = requests.get(
            f"{_AF_API}/{uniprot_id.upper()}", timeout=30
        )
        if api_resp.ok and api_resp.json():
            entry = api_resp.json()[0]
            pdb_url: str = entry["pdbUrl"]
            if fragment != 1:
                pdb_url = pdb_url.replace("-F1-", f"-F{fragment}-")
            filename = pdb_url.rsplit("/", 1)[-1]
        else:
            log.warning(
                "API returned no entry for %s (status %s); "
                "probing file versions %s directly.",
                uniprot_id, api_resp.status_code, _FALLBACK_VERSIONS,
            )
            for v in _FALLBACK_VERSIONS:
                candidate = (
                    f"AF-{uniprot_id.upper()}-F{fragment}-model_v{v}.pdb"
                )
                probe = requests.head(
                    f"{_AF_BASE}/{candidate}", timeout=20
                )
                if probe.ok:
                    filename, pdb_url = candidate, f"{_AF_BASE}/{candidate}"
                    log.info("Found %s at version v%d.", uniprot_id, v)
                    break
            else:
                raise FileNotFoundError(
                    f"No AlphaFold-DB file found for {uniprot_id!r} "
                    f"(tried versions {_FALLBACK_VERSIONS})."
                )
    else:
        filename = f"AF-{uniprot_id.upper()}-F{fragment}-model_v{version}.pdb"
        pdb_url = f"{_AF_BASE}/{filename}"

    dest = out_dir / filename
    log.info("Downloading %s", pdb_url)
    return _resume_get(pdb_url, dest)


def download_proteome(species: str, out_dir: Path) -> Path:
    """
    Fetch the full AlphaFold-DB v4 proteome tarball for *species*.

    Parameters
    ----------
    species:
        Case-insensitive species code.  Supported values:

        * ``"HUMAN"`` — H. sapiens (UP000005640, ~11 GB)
        * ``"MOUSE"`` — M. musculus (UP000000589)
        * ``"RAT"``   — R. norvegicus (UP000002494)
        * ``"YEAST"`` — S. cerevisiae (UP000002311)

    out_dir:
        Directory where the tarball will be saved.

    Returns
    -------
    Path
        Path to the ``.tar`` file; download is resumed if interrupted.

    Notes
    -----
    URL pattern::

        https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
            UP000005640_9606_HUMAN_v4.tar

    Do **not** call this function during smoke tests — the H. sapiens
    tarball is approximately 11 GB compressed.
    """
    species = species.upper()
    if species not in _PROTEOME_IDS:
        raise ValueError(
            f"Unknown species {species!r}. "
            f"Supported: {sorted(_PROTEOME_IDS)}"
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proteome_id = _PROTEOME_IDS[species]
    filename = f"{proteome_id}_v4.tar"
    dest = out_dir / filename

    url = f"{_PROTEOME_BASE}/{filename}"
    log.info("Downloading proteome tarball %s", url)
    return _resume_get(url, dest)


# ---------------------------------------------------------------------------
# Minimal inline PDB parser for the smoke test
# (parse.py, which uses gemmi, is implemented in Days 3-4)
# ---------------------------------------------------------------------------


def _quick_pdb_stats(pdb_path: Path) -> tuple[int, float]:
    """
    Count Cα residues and compute mean B-factor (= pLDDT) from a PDB file.

    Reads only ``ATOM`` records with atom name ``CA``; one unique
    (chain, resSeq, iCode) tuple per residue.  No external dependencies.

    Returns (n_residues, mean_pLDDT).
    """
    seen: set[tuple[str, int, str]] = set()
    bfactors: list[float] = []

    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            try:
                bfactors.append(float(line[60:66]))
            except ValueError:
                pass

    n = len(bfactors)
    mean_b = sum(bfactors) / n if n else float("nan")
    return n, mean_b


# ---------------------------------------------------------------------------
# Entry point — smoke test on the three prototype proteins
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    root = Path(__file__).resolve().parent.parent
    prototypes_dir = root / "data" / "prototypes"

    PROTOTYPES: list[tuple[str, str]] = [
        ("P15121",     "Ordered    (H. sapiens)  "),
        ("A0A0G2L439", "Disordered (Zebrafish)   "),
        ("Q9VQS4",     "Mixed      (D. melano.)  "),
    ]

    print("\nDownloading three prototype proteins from AlphaFold-DB v4 …\n")
    rows: list[tuple[str, str, int, float, float]] = []

    for uid, label in PROTOTYPES:
        path = download_structure(uid, prototypes_dir)
        n_res, mean_plddt = _quick_pdb_stats(path)
        size_kb = path.stat().st_size / 1024
        rows.append((label, uid, n_res, mean_plddt, size_kb))

    # ── print table ────────────────────────────────────────────────────────
    col = dict(type=28, uid=12, n=5, plddt=10, size=10)
    header = (
        f"{'Type':<{col['type']}}  {'UniProt':<{col['uid']}}  "
        f"{'N_aa':>{col['n']}}  {'mean pLDDT':>{col['plddt']}}  "
        f"{'Size (KB)':>{col['size']}}  Check"
    )
    sep = "-" * len(header)
    print(header)
    print(sep)

    failed: list[str] = []
    for label, uid, n, mplddt, kb in rows:
        ok = 0.0 <= mplddt <= 100.0
        flag = "OK" if ok else "WARN — out of range!"
        if not ok:
            failed.append(uid)
        print(
            f"{label:<{col['type']}}  {uid:<{col['uid']}}  "
            f"{n:>{col['n']}}  {mplddt:>{col['plddt']}.2f}  "
            f"{kb:>{col['size']}.1f}  {flag}"
        )

    print()
    print("Convention check: pLDDT should be in [0, 100] (AlphaFold-DB B-factor column).")

    if failed:
        print(f"FAILED: pLDDT out of range for {failed}")
        sys.exit(1)

    print("All checks passed. Smoke test OK.")
