# Q4 — pLDDT Fragmentation of AlphaFold Predictions at Proteome Scale

Course project for **Algorithms in Structural Bioinformatics** (Academic Year 2025–2026).

This repository reproduces **Q4** of Cazals & Sarti (2025) — pLDDT-based fragmentation analysis of AlphaFold reconstructions — at the scale of the complete *H. sapiens* proteome (23,391 AlphaFold-DB fragment files), and extends it with a **Q1 arity cross-correlation** analysis.

**What this measures.** AlphaFold tags every residue with a confidence score (pLDDT, 0–100) — high in well-folded regions, low in disordered or flexible ones. Sweeping the residue chain in order of *decreasing* pLDDT and tracking how it breaks into connected high-confidence components (a path-graph filtration solved with Union-Find under the Elder Rule) yields a persistence diagram. From it we derive the **persistence entropy** *H*<sub>p</sub> (how fragmented the confidence profile is) and *f*⁺<sub>cp</sub> (the fraction of long-lived components); a second-layer count of **persistent local maxima** (PLM) of the connected-component curve *N*<sub>cc</sub>(pLDDT) counts the distinct confident stretches. Together these flag multi-domain and partially-disordered proteins — confident segments punctuated by low-confidence linkers — purely from the 1-D pLDDT signal.

The path-graph filtration and persistence diagram (Union-Find with the Elder Rule) are implemented **from scratch** in Python. The only TDA library used is **GUDHI**, and only for the second-layer Persistent Local Maxima (PLM) computation on the *N*<sub>cc</sub> curve.

**Supervisor:** Prof. I. Emiris · **Co-advisor:** P. Rigas

---

## Paper references

F. Cazals & E. Sarti (2025).
*AlphaFold predictions on whole genomes at a glance: a coherent view on packing properties, pLDDT values, and disordered regions.*
bioRxiv 2024.11.16.623929 **v5** (posted April 10, 2025) — primary reference.
[https://www.biorxiv.org/content/10.1101/2024.11.16.623929v5](https://www.biorxiv.org/content/10.1101/2024.11.16.623929v5)

F. Cazals & E. Sarti (2025).
*AlphaFold predictions on whole genomes at a glance.*
bioRxiv 2024.11.16.623929 **v4** (posted January 3, 2025) — earlier version, consulted for the null model and Example 1.
[https://www.biorxiv.org/content/10.1101/2024.11.16.623929v4](https://www.biorxiv.org/content/10.1101/2024.11.16.623929v4)

---

## Results at a glance

H. sapiens proteome (AlphaFold-DB v4, 23,391 fragments):

| Metric | Paper | Batched default | Sequential (integer) | Raw float (sequential) |
|---|---|---|---|---|
| Pearson *r*(*f*⁺<sub>cp</sub>, *H*<sub>p</sub>) | ≈ 0.97 | **0.849** | 0.849 | 0.864 |
| Fig. 8 candidates (*n* ≥ 200, *H*<sub>p</sub> ≥ 0.25, PLM ≥ 3 at *t*<sub>p</sub> = 0.025) | ≈ 86 | **43** | 164 | 183 |
| Q1×Q4 arity centroid of candidate set | — | **(8, 15)** | (7, 15) | (7, 15) |
| Q1×Q4 enrichment at centroid bin | — | **10.88×** (*p* = 2.64×10⁻³) | 5.92× (*p* = 1.61×10⁻³) | 5.92× |

Prototype proteins (Figure 3 of the paper, integer pLDDT in all three columns of our pipeline):

| Protein | Class | Paper *H*<sub>p</sub> | Our *H*<sub>p</sub> | Δ |
|---|---|---|---|---|
| P15121 | Ordered | ≈ 0.04 | 0.147 | +0.11 |
| A0A0G2L439 | Disordered | ≈ 0.27 | 0.433 | +0.16 |
| Q9VQS4 | Mixed | ≈ 0.28 | 0.425 | +0.15 |

The residual gap to the paper's per-protein *H*<sub>p</sub> and Fig. 8 numbers is attributed to **two non-exclusive causes** that we cannot separate without the authors' original inputs: (1) **dataset drift** between the AlphaFold-DB v4 snapshot the paper authors used in 2024 and the v4 tarball EBI currently distributes (v4 has been periodically rebuilt as AlphaFold 2.3 received sequence-database refreshes); and (2) the paper's **undisclosed pLDDT pre-processing** — no single rounding width simultaneously reproduces the prototype *H*<sub>p</sub> values and the Fig. 8 candidate count, so we present integer pLDDT as our best-justified default rather than a claim to have recovered the paper's exact convention. With neither the authors' frozen snapshot nor their pre-processing code available, the two contributions cannot be disentangled. See `report/main.pdf` Section VI for the full discussion (including a side-by-side *N*<sub>cc</sub> curve overlay, the proteome-wide PLM histogram, and the Fig. 8 candidate count bar chart). The paper's 86-candidate count sits squarely between our batched-default 43 and our sequential-ablation 164 (geometric mean √(43·164) ≈ 84).

### Cross-organism generalisation

Re-running the **identical** integer-pLDDT, batched pipeline on **twelve** further reference proteomes (AlphaFold-DB v6) spanning **all three domains of life** — chosen to break the size–clade confound by pairing large non-mammalian eukaryotes against small proteomes from every clade — gives the *f*⁺<sub>cp</sub>–*H*<sub>p</sub> correlation across 13 organisms (eukaryotes first, prokaryotes below):

| Organism | Fragments | Pearson *r* | median *H*<sub>p</sub> | Fig. 8 candidates |
|---|---|---|---|---|
| *H. sapiens* | 23,586 | 0.850 | 0.317 | 43 |
| *M. musculus* | 21,452 | 0.862 | 0.313 | 37 |
| *R. norvegicus* | 22,152 | 0.854 | 0.321 | 38 |
| *D. rerio* (vertebrate) | 26,290 | 0.876 | 0.314 | 66 |
| *D. melanogaster* (invertebrate) | 13,461 | 0.822 | 0.320 | 26 |
| *C. elegans* (invertebrate) | 19,700 | 0.846 | 0.322 | 42 |
| *A. thaliana* (plant) | 27,402 | 0.864 | 0.316 | 65 |
| *S. cerevisiae* (fungus) | 6,055 | 0.816 | 0.311 | 11 |
| *P. falciparum* (protist) | 5,168 | 0.764 | 0.332 | 28 |
| *P. aeruginosa* (bacterium) | 5,555 | 0.798 | 0.283 | 4 |
| *E. coli* (bacterium) | 4,370 | 0.745 | 0.286 | 1 |
| *M. tuberculosis* (bacterium) | 3,991 | 0.842 | 0.288 | 4 |
| *M. jannaschii* (archaeon) | 1,773 | 0.786 | 0.283 | 2 |

(All v6; the *H. sapiens* row here is the v6 counterpart of the *r* = 0.849 v4 figure above.) The evidence reads in **two layers**. The fragmentation **phenomenon** — the *level* of the statistics — is essentially identical across all **nine eukaryotes** (median *f*⁺<sub>cp</sub> 0.194 ± 0.007, median *H*<sub>p</sub> 0.318 ± 0.006; the protist is the most disorder-rich of all); the **four prokaryotes** are structurally flatter (median *f*⁺<sub>cp</sub> 0.153 ± 0.004, *H*<sub>p</sub> 0.285 ± 0.002) but still show the coupling, so the signal is universal across the eukaryotic domain and attenuated, never absent, beyond it. The **strength** of the correlation (*r* = 0.745–0.876) is **governed by proteome size**: with the size–clade confound broken, per-organism *r* rises with fragment count at **Spearman 0.86 (*p* = 0.001, n = 13)** independently of clade — the three large non-mammalian eukaryotes (zebrafish, *Arabidopsis*, *C. elegans*) average *r* = 0.862, identical to the three mammals' 0.856. **Intrinsic-disorder content is excluded**: *r* is uncorrelated with median *H*<sub>p</sub> (Spearman 0.31, *p* = 0.31) and the most disordered proteome (*P. falciparum*) has the second-lowest *r* of all 13. *M. tuberculosis* (small but *r* = 0.84) is the one outlier; dropping it lifts the size association to Spearman 0.94. The arity cross-correlation now replicates **independently significantly** in zebrafish (11.2×, *p* = 2.5×10⁻³) and *C. elegans* (9.7×, *p* = 1.8×10⁻²); the four prokaryotes yield too few candidates (1–4) to test. See `report/main.pdf` Section VI and `results/figures/cross_organism_r_vs_size.png`.

---

## Pipeline defaults

Two algorithmic choices distinguish our default pipeline from the most literal reading of the paper, both of them needed to reproduce the paper's null-model behaviour and qualitative findings:

| Choice | Default | Flag | Rationale |
|---|---|---|---|
| **pLDDT discretisation** | round to integers in [0, 100] | `discretise=True` | AlphaFold itself reports per-residue confidence as an integer. The raw-float baseline matches the paper's null *H*<sub>p</sub> at *n* = 1,000 only by coincidence and drifts to ≈ 0.41 at *n* = 10,000 (paper: 0.47); integer pLDDT reproduces all three null *H*<sub>p</sub> values across the paper's three sample sizes. |
| **Batched insertion** | residues at the same integer pLDDT level inserted together | `batched=True` | Matches the paper's statement that "pLDDT values come in batches". The persistence diagram is bit-for-bit identical to the sequential one-residue-per-step variant under the Elder Rule, but the *N*<sub>cc</sub> curve is piecewise constant within each plateau, which significantly reduces PLM counts and the Fig. 8 candidate count. |

Both flags can be set to `False` on `build_pd_and_ncc(plddt, discretise=..., batched=...)` to access the raw-float and sequential baselines as ablations.

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
  run_organism.py      – full pipeline for any model organism (13 supported, all 3 domains of life)
  compare_organisms.py – cross-organism comparison table + Pearson-r bar chart + r-vs-size scatter
  plots.py             – figure7 / figure8_scatter / figure8_panels helpers
tests/                 – pytest unit tests for every algorithm module (57 tests)
notebooks/             – seven Jupyter notebooks (walkthrough → proteome figures,
                         incl. per-organism copies of notebook 04)
report/                – LaTeX source (IEEEtran, gitignored) + compiled PDF (tracked)
data/
  prototypes/          – prototype PDB files for the 3 proteins (bundled, ~1.2 MB)
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

All 57 tests should pass. The suite covers `filtration`, `statistics`, `null_model`, `plm`, `parse`, and `arity`, including regression tests for the three prototype proteins under the integer + batched default and for the inverse-CDF arity signature (paper Def. 2).

---

## Reproducing the results

### Step 1 — Prototype proteins (bundled)

The PDB files for the three prototype proteins are committed under
`data/prototypes/` (~1.2 MB), so prototype validation (Step 2) needs no
download. `python src/download.py` is kept as an optional refresh, but note
that AlphaFold-DB has since removed some individual files (A0A0G2L439 is no
longer available; P15121/Q9VQS4 now resolve to v6), so it may 404.

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

### Step 3 — Full proteome run (H. sapiens v4)

On first run, `run_full_proteome.py` downloads the *H. sapiens* AlphaFold-DB
v4 tarball (`UP000005640_9606_HUMAN_v4.tar`, ~5 GB) into `data/hsapiens/`
(resumable if interrupted), then streams it. The two follow-up scripts reuse
that same local tarball, so run them in order:

```bash
python src/run_full_proteome.py    # writes data/results/proteome_full.parquet
python src/run_arity.py            # writes data/results/proteome_full_arity.parquet
python src/run_tp_ablation.py      # writes data/results/proteome_full_tp_ablation.parquet
```

After the one-time download, each pass takes about 30 s on 6 worker processes.

### Step 4 — v6 replication and cross-organism extension

These reproduce the generalisation results in the discussion (Section VI of
`report/main.pdf`). Each runner auto-downloads its AlphaFold-DB **v6** tarball
into `data/<organism>/` (resumable) and writes the parquets that
`compare_organisms.py` consumes:

```bash
# H. sapiens v6 (~5 GB tarball): main stats, then the matching arity pass
python src/run_v6_proteome.py
python -c "import sys; sys.path.insert(0, 'src'); import run_arity; \
  run_arity.run(tar_path='data/hsapiens/UP000005640_9606_HUMAN_v6.tar', \
                out_path='data/results/proteome_v6_arity.parquet')"

# Mouse / rat / yeast v6 (each writes full + tp_ablation + arity parquets)
python src/run_organism.py MOUSE --version v6
python src/run_organism.py RAT   --version v6
python src/run_organism.py YEAST --version v6

# P. aeruginosa — bacterial out-of-clade probe (~640 MB v6 tarball, ~5.5k proteins)
python src/run_organism.py PSEAE --version v6

# D. melanogaster — invertebrate (non-mammalian animal) probe (~2.3 GB v6 tarball, ~13.5k proteins)
python src/run_organism.py DROME --version v6

# Size×clade confound-breakers spanning all three domains of life. Tarballs are
# large; with limited disk, delete each after its run (parquets are <1 MB):
#   for SP in DANRE ARATH CAEEL PLAF7 ECOLI MYCTU METJA; do
#     python src/run_organism.py $SP --version v6 && rm data/${SP:l}/*.tar
#   done
python src/run_organism.py DANRE --version v6   # zebrafish — vertebrate (~5.0 GB)
python src/run_organism.py ARATH --version v6   # A. thaliana — plant (~3.9 GB)
python src/run_organism.py CAEEL --version v6   # C. elegans — invertebrate (~2.8 GB)
python src/run_organism.py PLAF7 --version v6   # P. falciparum — protist (~1.2 GB)
python src/run_organism.py ECOLI --version v6   # E. coli — bacterium (~0.5 GB)
python src/run_organism.py MYCTU --version v6   # M. tuberculosis — bacterium (~0.4 GB)
python src/run_organism.py METJA --version v6   # M. jannaschii — archaeon (~0.2 GB)

# Cross-organism comparison table + figures (bar chart + r-vs-log10(N) scatter)
python src/compare_organisms.py    # writes results/figures/cross_organism_{pearson,r_vs_size}.png
```

The explicit `python -c` for the v6 arity pass is needed because `run_arity.py`
defaults to the v4 tarball and has no command-line override; the helper calls
its `run()` entry point directly. `compare_organisms.py` uses whichever
organism parquets are present, so you can run a subset.

### Notebooks

For a guided walkthrough open the notebooks in order:

```bash
jupyter notebook notebooks/
```


| Notebook                              | Content                                       |
| ------------------------------------- | --------------------------------------------- |
| `01_single_protein_walkthrough.ipynb` | Filtration step-by-step on one protein        |
| `02_null_model.ipynb`                 | Random pLDDT baseline and Conjecture 1        |
| `03_prototypes_figure3.ipynb`         | Reproduce Figure 3 (three prototype proteins) |
| `04_proteome_figures7_8.ipynb`        | Reproduce Figures 7 and 8 at proteome scale   |
| `04_proteome_figures7_8_mouse.ipynb`  | Figures 7/8 for *M. musculus* (v6)            |
| `04_proteome_figures7_8_rat.ipynb`    | Figures 7/8 for *R. norvegicus* (v6)          |
| `04_proteome_figures7_8_yeast.ipynb`  | Figures 7/8 for *S. cerevisiae* (v6)          |


---

## Dependencies

All runtime and test dependencies are listed in `requirements.txt`.
Key packages: `numpy`, `scipy`, `gemmi`, `pandas`, `pyarrow`, `matplotlib`,
`gudhi`, `requests`, `tqdm`, `pytest`.

---

*Giorgos Boulogeorgos — Algorithms in Structural Bioinformatics, 2025–2026*
