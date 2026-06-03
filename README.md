# Q4 — pLDDT Fragmentation of AlphaFold Predictions at Proteome Scale

Course project for **Algorithms in Structural Bioinformatics** (Academic Year 2025–2026).

This repository reproduces **Q4** of Cazals & Sarti (2025) — pLDDT-based fragmentation analysis of AlphaFold reconstructions — at the scale of the complete *H. sapiens* proteome (23,391 AlphaFold-DB fragment files), and extends it with a **Q1 arity cross-correlation** analysis.

The path-graph filtration and persistence diagram (Union-Find with the Elder Rule) are implemented **from scratch** in Python. The only TDA library used is **GUDHI**, and only for the second-layer Persistent Local Maxima (PLM) computation on the *N*<sub>cc</sub> curve.

**Supervisor:** Prof. I. Emiris · **Co-advisor:** P. Rigas

---

## Paper reference

F. Cazals & A. Sarti (2025).
*AlphaFold predictions on whole genomes at a glance: a coherent view on packing properties, pLDDT values, and disordered regions.*
bioRxiv 2024.11.16.623929 v4.
[https://doi.org/10.1101/2024.11.16.623929](https://doi.org/10.1101/2024.11.16.623929)

---

## Pipeline defaults

Two algorithmic choices distinguish our default pipeline from the most literal reading of the paper, both of them needed to reproduce the paper's null-model behaviour and qualitative findings:

| Choice | Default | Flag | Rationale |
|---|---|---|---|
| **pLDDT discretisation** | round to integers in [0, 100] | `discretise=True` | AlphaFold itself reports per-residue confidence as an integer. The raw-float baseline matches the paper's null *H*<sub>p</sub> at *n* = 1,000 only by coincidence and drifts to ≈ 0.41 at *n* = 10,000 (paper: 0.47); integer pLDDT reproduces all three null *H*<sub>p</sub> values across the paper's three sample sizes. |
| **Batched insertion** | residues at the same integer pLDDT level inserted together | `batched=True` | Matches the paper's statement that "pLDDT values come in batches". The persistence diagram is bit-for-bit identical to the sequential one-residue-per-step variant under the Elder Rule, but the *N*<sub>cc</sub> curve is piecewise constant within each plateau, which significantly reduces PLM counts and the Fig. 8 candidate count. |

Both flags can be set to `False` on `build_pd_and_ncc(plddt, discretise=..., batched=...)` to access the raw-float and sequential baselines as ablations.

---

## Results at a glance

H. sapiens proteome (AlphaFold-DB v4, 23,391 fragments):

| Metric | Paper | Batched default | Sequential (integer) | Raw float (sequential) |
|---|---|---|---|---|
| Pearson *r*(*f*⁺<sub>cp</sub>, *H*<sub>p</sub>) | ≈ 0.97 | **0.849** | 0.849 | 0.864 |
| Fig. 8 candidates (*n* ≥ 200, *H*<sub>p</sub> ≥ 0.25, PLM ≥ 3 at *t*<sub>p</sub> = 0.025) | ≈ 86 | **43** | 164 | 183 |
| Q1×Q4 arity centroid of candidate set | — | **(8, 15)** | (7, 15) | (7, 15) |
| Q1×Q4 enrichment at centroid bin | — | **10.88×** (*p* = 2.64×10⁻³) | 5.92× (*p* = 1.61×10⁻³) | 5.92× |

Prototype proteins (Figure 1 of the paper, integer pLDDT in all three columns of our pipeline):

| Protein | Class | Paper *H*<sub>p</sub> | Our *H*<sub>p</sub> | Δ |
|---|---|---|---|---|
| P15121 | Ordered | ≈ 0.04 | 0.147 | +0.11 |
| A0A0G2L439 | Disordered | ≈ 0.27 | 0.433 | +0.16 |
| Q9VQS4 | Mixed | ≈ 0.28 | 0.425 | +0.15 |

The residual gap to the paper's per-protein and Fig. 8 numbers is attributed to dataset drift between the AlphaFold-DB v4 snapshot used by the paper authors and the v4 tarball currently distributed by EBI; see `report/main.pdf` Section VI for the full discussion (including a side-by-side *N*<sub>cc</sub> curve overlay, the proteome-wide PLM histogram, and the Fig. 8 candidate count bar chart). The paper's 86-candidate count sits squarely between our batched-default 43 and our sequential-ablation 164 (geometric mean √(43·164) ≈ 84).

---

## Repository layout

```
src/
  parse.py             – extract per-residue pLDDT from AlphaFold PDB files (gemmi)
  filtration.py        – build_pd_and_ncc(plddt, discretise=True, batched=True)
                          path-graph filtration + Elder Rule persistence diagram
  statistics.py        – f_cp_plus, mean_persistence, persistence_entropy, ncc_max
  plm.py               – second-layer PLM via GUDHI CubicalComplex
  arity.py             – Cα-packing arity signatures (a25, a75)
  download.py          – fetch AlphaFold-DB structures
  process_proteome.py  – random subset pipeline (default N=500)
  run_full_proteome.py – streaming pipeline over the v4 tarball
  run_v6_proteome.py   – same pipeline against the v6 tarball
  run_tp_ablation.py   – three-threshold PLM ablation at t_p ∈ {0.020, 0.025, 0.030}
  run_arity.py         – arity signatures proteome-wide
  plots.py             – figure7 / figure8_scatter / figure8_panels helpers
tests/                 – pytest unit tests for every algorithm module (56 tests)
notebooks/             – four Jupyter notebooks (walkthrough → proteome figures)
report/                – LaTeX source (IEEEtran, gitignored) + compiled PDF (tracked)
slides/                – Beamer LaTeX source + PDF + PPTX export
data/
  prototypes/          – three prototype PDB files (downloaded by download.py)
  hsapiens/            – AlphaFold-DB proteome tarballs (~5 GB each, not committed)
  results/             – parquet files produced by the pipeline
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

All 56 tests should pass. The suite covers `filtration`, `statistics`, `null_model`, `plm`, `parse`, and `arity`, including regression tests for the three prototype proteins under the integer + batched default.

---

## Reproducing the results

### Step 1 — Download the three prototype proteins

```bash
python src/download.py
```

Downloads ~500 KB of PDB files into `data/prototypes/`. No proteome data needed for prototype validation.

### Step 2 — Validate prototypes

```bash
python -c "
import sys; sys.path.insert(0, 'src')
from parse import parse_alphafold_pdb
from filtration import build_pd_and_ncc
from statistics import persistence_entropy, f_cp_plus

PROTOS = [
    ('P15121',     0.04, 'data/prototypes/AF-P15121-F1-model_v4.pdb'),
    ('A0A0G2L439', 0.27, 'data/prototypes/AF-A0A0G2L439-F1-model_v4.pdb'),
    ('Q9VQS4',     0.28, 'data/prototypes/AF-Q9VQS4-F1-model_v6.pdb'),
]
for uid, paper_hp, path in PROTOS:
    plddt = parse_alphafold_pdb(path)['plddt']
    pd, _ = build_pd_and_ncc(plddt)            # integer + batched defaults
    print(f'{uid:12s}  H_p={persistence_entropy(pd):.3f}  '
          f'f+_cp={f_cp_plus(pd):.3f}  (paper H_p ~ {paper_hp})')
"
```

Expected output under the integer + batched default: `H_p ≈ 0.147 / 0.433 / 0.425`. See *Results at a glance* above for the explanation of the residual gap to the paper.

### Step 3 — Full proteome run (requires ~5 GB download)

Download the *H. sapiens* AlphaFold-DB v4 tarball from
[https://alphafold.ebi.ac.uk/download](https://alphafold.ebi.ac.uk/download) into `data/hsapiens/`, then:

```bash
python src/run_full_proteome.py    # writes data/results/proteome_full.parquet
python src/run_arity.py            # writes data/results/proteome_full_arity.parquet
python src/run_tp_ablation.py      # writes data/results/proteome_full_tp_ablation.parquet
```

Each run takes about 30 s on 8 worker processes against a local tarball.

### Notebooks

For a guided walkthrough open the notebooks in order:

```bash
jupyter notebook notebooks/
```


| Notebook                              | Content                                       |
| ------------------------------------- | --------------------------------------------- |
| `01_single_protein_walkthrough.ipynb` | Filtration step-by-step on one protein        |
| `02_null_model.ipynb`                 | Random pLDDT baseline and Conjecture 1        |
| `03_prototypes_figure3.ipynb`         | Reproduce Figure 1 (three prototype proteins) |
| `04_proteome_figures7_8.ipynb`        | Reproduce Figures 7 and 8 at proteome scale   |


---

## Dependencies

All runtime and test dependencies are listed in `requirements.txt`.
Key packages: `numpy`, `scipy`, `gemmi`, `pandas`, `pyarrow`, `matplotlib`,
`gudhi`, `requests`, `tqdm`, `pytest`.

---

*Giorgos Boulogeorgos — Algorithms in Structural Bioinformatics, 2025–2026*
