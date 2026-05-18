# Q4 — pLDDT Fragmentation of AlphaFold Predictions at Proteome Scale

Reproduction of **Q4** of Cazals & Sarti (2025) applied to the complete
H. sapiens proteome (~20 000 AlphaFold-DB v4 predictions), plus a Q1
arity-map cross-correlation extension.

## Paper reference

F. Cazals & A. Sarti (2025). *AlphaFold predictions on whole genomes at a
glance: a coherent view on packing properties, pLDDT values, and disordered
regions*. bioRxiv 2024.11.16.623929 v4.
<https://doi.org/10.1101/2024.11.16.623929>

## Quickstart

```bash
# 1. Clone
git clone <repo-url> project_q4 && cd project_q4

# 2. Create environment (conda recommended — gudhi is on conda-forge)
conda create -n q4 python=3.11
conda activate q4
conda install -c conda-forge gudhi
pip install -r requirements.txt

# 3. Download the three prototype proteins and run the smoke test
python src/download.py
```

The smoke test downloads three small PDB files (~500 KB each) and prints a
table confirming residue counts and that pLDDT values (B-factor column) lie
in [0, 100]. It does **not** download the 11 GB proteome tarball.

## Repository layout

`src/` holds the analysis pipeline as importable modules: `download.py`
fetches AlphaFold-DB structures; `parse.py` extracts per-residue pLDDT from
PDB files; `filtration.py` implements the path-graph filtration and persistence
diagram via Union-Find with the Elder Rule; `statistics.py` computes the five
summary statistics (f⁺_cp, p̄, H_p, N_cc^max, PLM); `null_model.py` simulates
the random-pLDDT baseline (Example 1, Conjecture 1 of the paper); `plm.py`
wraps gudhi's CubicalComplex for the second-layer PLM persistence;
`process_proteome.py` runs the batch pipeline over the full proteome and writes
`results/proteome_stats.parquet`; `plots.py` generates Figures 3, 7, and 8.
Exploratory analyses live in `notebooks/`. The LaTeX report source is in
`report/` and presentation decks are in `slides/`.
