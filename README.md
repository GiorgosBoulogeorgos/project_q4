# Q4 — pLDDT Fragmentation of AlphaFold Predictions at Proteome Scale

![Python](https://img.shields.io/badge/python-3.11%2B-blue)

Course project for **Algorithms in Structural Bioinformatics** (Academic Year 2025–2026).

This repository reproduces **Q4** of Cazals & Sarti (2025) — pLDDT-based fragmentation analysis of AlphaFold reconstructions — at the scale of the complete *H. sapiens* proteome (~23 000 AlphaFold-DB fragment files), and extends it with a **Q1 arity cross-correlation** analysis.

The core algorithm (path-graph filtration + persistence diagram via Union-Find with the Elder Rule) is implemented **from scratch** in Python, without any TDA library, following the paper exactly.

**Supervisor:** Prof. I. Emiris &nbsp;·&nbsp; **Co-advisor:** P. Rigas

---

## Paper reference

F. Cazals & A. Sarti (2025).
*AlphaFold predictions on whole genomes at a glance: a coherent view on packing properties, pLDDT values, and disordered regions.*
bioRxiv 2024.11.16.623929 v4.
<https://doi.org/10.1101/2024.11.16.623929>

---

## Results at a glance

| Target (paper) | This reproduction |
|---|---|
| Pearson r(f⁺_cp, H_p) ≈ 0.97 | r = 0.864 (gap explained by float-precision pLDDT in v4) |
| ~86 fragmentation candidates at t_p = 0.025 | 183 candidates (same precision-shift explanation) |
| H_p ≈ 0.04 / 0.27 / 0.28 for ordered / disordered / mixed prototypes | Reproduced within ±0.03 |

---

## Repository layout

```
src/
  parse.py            – extract per-residue pLDDT from AlphaFold PDB files
  filtration.py       – path-graph filtration + Elder Rule persistence diagram
  statistics.py       – five summary statistics (f⁺_cp, p̄, H_p, N_cc_max, PLM)
  null_model.py       – random-pLDDT baseline (Example 1 / Conjecture 1)
  plm.py              – second-layer PLM via gudhi CubicalComplex
  arity.py            – Cα-packing arity signatures
  download.py         – fetch AlphaFold-DB structures
  process_proteome.py – batch pipeline → results/proteome_stats.parquet
  run_full_proteome.py / run_v6_proteome.py / run_arity.py / run_tp_ablation.py
  plots.py            – reproduce Figures 3, 7, 8 from the paper
tests/                – pytest unit tests for every algorithm module
notebooks/            – four Jupyter notebooks (walkthrough → proteome figures)
report/               – LaTeX source (IEEEtran) + compiled PDF
slides/               – Beamer LaTeX source + compiled PDF + PPTX export
data/
  prototypes/         – three prototype PDB files (downloaded by download.py)
  hsapiens/           – AlphaFold-DB proteome tarballs (not committed, ~11 GB)
results/              – parquet files produced by the pipeline
```

---

## Setup

### Prerequisites

- Python 3.11+
- [conda](https://docs.conda.io/) (recommended — `gudhi` is easiest via conda-forge)

### 1. Clone

```bash
git clone https://github.com/GiorgosBoulogeorgos/project_q4.git
cd project_q4
```

### 2. Create environment

**With conda (recommended):**

```bash
conda create -n q4 python=3.11
conda activate q4
conda install -c conda-forge gudhi
pip install -r requirements.txt
```

**With venv only** (if `pip install gudhi` works on your platform):

```bash
python3.11 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### 3. Verify the install

```bash
pytest tests/ -q
```

All tests should pass. The test suite covers the filtration, statistics,
null model, PLM, parse, and arity modules, including sanity checks against
the three prototype proteins.

---

## Reproducing the results

### Step 1 — Download the three prototype proteins

```bash
python src/download.py
```

Downloads ~500 KB of PDB files into `data/prototypes/`. No proteome data
needed for the prototype validation.

### Step 2 — Validate prototypes (Figure 3)

```bash
python -c "
from src.parse import parse_plddt
from src.filtration import compute_persistence_diagram
from src.statistics import compute_all_stats
import pathlib

for uid, expected_hp in [('P15121', 0.04), ('A0A0G2L439', 0.27), ('Q9VQS4', 0.28)]:
    pdbs = sorted(pathlib.Path('data/prototypes').glob(f'*{uid}*'))
    pds = [compute_persistence_diagram(parse_plddt(p)) for p in pdbs]
    stats = compute_all_stats(pds)
    print(f'{uid}  H_p={stats[\"H_p\"]:.3f}  (expected ~{expected_hp})')
"
```

### Step 3 — Full proteome run (requires ~11 GB download)

Download the *H. sapiens* AlphaFold-DB v4 tarball from
<https://alphafold.ebi.ac.uk/download> into `data/hsapiens/`, then:

```bash
python src/run_full_proteome.py          # writes results/proteome_stats.parquet
python src/run_arity.py                  # writes results/arity_stats.parquet
python src/plots.py                      # regenerates figures 3, 7, 8
```

### Notebooks

For a guided walkthrough open the notebooks in order:

```bash
jupyter notebook notebooks/
```

| Notebook | Content |
|---|---|
| `01_single_protein_walkthrough.ipynb` | Filtration step-by-step on one protein |
| `02_null_model.ipynb` | Random pLDDT baseline and Conjecture 1 |
| `03_prototypes_figure3.ipynb` | Reproduce Figure 3 (three prototype proteins) |
| `04_proteome_figures7_8.ipynb` | Reproduce Figures 7 and 8 at proteome scale |

---

## Dependencies

All runtime and test dependencies are listed in `requirements.txt`.
Key packages: `numpy`, `scipy`, `gemmi`, `pandas`, `pyarrow`, `matplotlib`,
`gudhi`, `requests`, `tqdm`, `pytest`.

---

*Giorgos Boulogeorgos — Algorithms in Structural Bioinformatics, 2025–2026*
