# Project 10 (Q4 variant) — Project Plan

**Course:** Algorithms in Structural Bioinformatics
**Reference paper:** Cazals & Sarti, *AlphaFold predictions on whole genomes at a glance: a coherent view on packing properties, pLDDT values, and disordered regions*, bioRxiv 2024.11.16.623929 (v4, April 2025)
**Question tackled:** Q4 — pLDDT values and fragmentation of AlphaFold reconstructions
**Scope:** Full reproduction at proteome scale (H. Sapiens, ~20k AlphaFold predictions)
**Timeline:** 3 weeks

---

## 1. Goal

Reproduce the Q4 analysis of Cazals & Sarti at the scale of the H. Sapiens proteome. Concretely:

1. For each AlphaFold prediction, build the filtration **G_pLDDT** of the polypeptide path graph and compute its persistence diagram.
2. Derive the five statistics defined in the paper: f⁺_cp, mean persistence p̄, normalized persistence entropy H_p, N_cc^max, and PLM(t_ν).
3. Reproduce the paper's null random model (Example 1, Conjecture 1) and use it as a baseline.
4. Reproduce the qualitative behaviour on the three prototypical proteins (Figure 3 of the new paper):
   - **Ordered:** H. Sapiens AF-P15121, 316 a.a., expected H_p ≈ 0.04
   - **Disordered:** Zebrafish AF-A0A0G2L439, 449 a.a., expected H_p ≈ 0.27
   - **Mixed/fragmented:** D. Melanogaster AF-Q9VQS4, 781 a.a., expected H_p ≈ 0.28
5. Run the pipeline on the entire H. Sapiens proteome and reproduce:
   - **Figure 7:** scatter plot of f⁺_cp vs. H_p, expected Pearson r ≈ 0.97
   - **Figure 8:** 3D scatter (#a.a., H_p, f⁺_cp); query the fragmented subset {H_p ≥ 0.25, #a.a. ≥ 200, PLM(t_p=0.025) ≥ 3}. Paper reports 86 such structures.
6. Inspect and visualize a sample of the fragmented predictions.

---

## 2. Conceptual summary of Q4 (for the report and slides)

### 2.1 The path-graph filtration G_pLDDT

A polypeptide chain of n residues is modelled as a path graph: n vertices, n−1 edges connecting consecutive residues. Each vertex carries a pLDDT score in [0, 100]. The filtration G_pLDDT inserts residues in order of **decreasing pLDDT** (equivalently: increasing −pLDDT). When residue *i* is inserted, the edges (*i*−1, *i*) and (*i*, *i*+1) are added if those neighbours are already present.

As insertions proceed, the number of connected components evolves: a function N_cc(pLDDT) defined on the sweep parameter.

### 2.2 Persistence diagram of G_pLDDT

Each component is born when its first residue is inserted (pLDDT = b) and dies when it merges with an older component (pLDDT = d, with d < b since pLDDT is decreasing). The pair (b, d) is a point in the persistence diagram (PD). Points with b = d are **accretions** (a residue inserted adjacent to a just-formed component); points with b > d have positive persistence.

### 2.3 Summary statistics

| Statistic | Definition | Intuition |
|---|---|---|
| f⁺_cp | m′/(n−1), with m′ = # critical points with positive persistence | Fraction of "real" merges (not trivial accretions). High → scrambled pLDDT along sequence. |
| p̄ | mean persistence Σ p(c_i) / m | Average lifetime of a component. |
| H_p | normalized Shannon entropy of persistence distribution | Diversity of lifetimes. Low → uniform; high → broad range, often fragmented protein. |
| N_cc^max | max of N_cc over the sweep | Peak component count. Null model conjecture: ≈ n/4. |
| PLM(t_ν) | # persistent local maxima of N_cc after Morse-Smale simplification at threshold t_ν = n · t_p | Real "peaks" in the N_cc curve, i.e. independent stretches that form before merging. |

### 2.4 Null model

A protein with random pLDDT ∼ U[0,100] gives, empirically:
- f⁺_cp ≈ 0.33 (independent of n)
- H_p ≈ 0.36 (n=100) → 0.47 (n=10,000)
- N_cc^max / n ≈ 1/4 (Conjecture 1)
- PD has a characteristic triangular envelope

Real proteins are *less* random than this: ordered proteins have far lower H_p and a near-empty PD; disordered proteins approach the null model; fragmented predictions sit somewhere in between with a wide spread of persistences.

---

## 3. Deliverables

1. **Introductory slides (PPTX)** — short deck (~8 slides) introducing the project, given before any work is done.
2. **Code repository** — Python scripts and a Jupyter notebook implementing the full pipeline.
3. **LaTeX report (PDF)** — methods, results, figures, discussion. Target: ~12-15 pages.
4. **Final presentation slides (PPTX)** — ~15 slides covering motivation, methods, results, discussion.

---

## 4. Technical stack

- **Python 3.11+**
- numpy, scipy — numerics
- gemmi — fast PDB/mmCIF parsing (proteome scale)
- pandas — tabular per-protein statistics
- matplotlib, seaborn — plots
- **gudhi** — second-layer persistence (for the PLM statistic only, applied to the 1D function N_cc)
- Path-graph filtration and Union-Find implemented from scratch
- tqdm — progress reporting
- requests, tarfile — proteome download

Justification for the hybrid implementation choice: the filtration of a path graph is most efficient as hand-rolled Union-Find with the Elder Rule; it is also more educational to implement directly. The second persistence computation, on the 1D function N_cc, has more edge cases (boundary handling, tie-breaking) and is cleanly handled by gudhi's `CubicalComplex`. This split is described explicitly in the methods section of the report.

---

## 5. Three-week timeline

### Week 1 — Foundations and single-protein analysis

**Days 1–2: Introductory slides + repo skeleton**
- Create the introductory PPTX (8 slides) for the professor.
- Set up the project repo: directory structure, `requirements.txt`, README.
- Write the data-download script for AlphaFold-DB (single-protein and full-proteome modes).

**Days 3–4: Core algorithms on a single protein**
- Implement `parse_alphafold_pdb()` extracting per-residue pLDDT (the B-factor column for AlphaFold predictions).
- Implement `build_filtration_and_pd()`: Union-Find with the Elder Rule, returning the PD as an array of (birth, death) pairs and the N_cc curve.
- Implement the five summary statistics.
- Implement the null random model (verify Conjecture 1 and the reported f⁺_cp ≈ 0.33).

**Day 5: Reproduce Figure 3**
- Download the three prototype proteins.
- Generate the full figure: 3D structure (or just a placeholder), N_cc(pLDDT) curve with cumulative residue fraction overlay, and persistence diagram.
- Sanity-check H_p values against the paper (0.04, 0.27, 0.28).

**Days 6–7: PLM and Morse-Smale simplification**
- Implement PLM(t_ν) via gudhi: feed the 1D N_cc array as a cubical complex, extract 0-dimensional persistence pairs, count pairs with persistence ≥ t_ν.
- Validate against the Bottom example of Figure 3 (D. Melanogaster, should yield PLM ≥ 2 at reasonable t_p).
- Write a small unit-test file covering all single-protein components.

**End of Week 1 milestone:** all algorithms work correctly on the three prototype proteins; introductory slides ready.

### Week 2 — Proteome-scale analysis

**Days 8–9: Proteome download and batch pipeline**
- Download the H. Sapiens AlphaFold-DB proteome tarball (~11 GB compressed; instruction in report on how to obtain it).
- Write `process_proteome.py` that streams through all structures, computes the five statistics per protein, and produces a single tidy `results.parquet` file.
- Add error handling for fragmented PDB entries (proteins with multiple "F1, F2, ..." fragments — decide whether to concatenate or treat separately; paper appears to treat them separately).

**Days 10–11: Reproduce Figure 7**
- Compute Pearson correlation between f⁺_cp and H_p across the proteome. Target: r ≈ 0.97.
- Plot the scatter with min-H_p and max-H_p insets showing the corresponding PDs and structures.
- Cross-check the cumulative pLDDT distribution against Table 1 of the paper (H. Sapiens row: 0.284 / 0.382 / 0.666 at pLDDT thresholds 50/70/90).

**Days 12–13: Reproduce Figure 8**
- Compute PLM for every protein at t_p = 0.025.
- 3D scatter (#a.a., H_p, f⁺_cp).
- Filter on (#a.a. ≥ 200, H_p ≥ 0.25, PLM ≥ 3). Count the surviving proteins — target: 86 (paper).
- Display three of them as in panels (B, C, D) of Figure 8.

**Day 14: Extension — cross-correlation with Q1 (arity map)**
- Implement the arity computation from Q1: for each Cα atom, count Cα neighbours within 10 Å.
- For each protein, compute the arity signature (arity values at the 25th and 75th percentiles of the per-residue arity distribution).
- Bin all proteins on the 2D arity map and overlay our Q4 statistics (mean H_p, mean PLM, fraction of fragmented proteins) per bin.
- Question: do the 86 "fragmented" proteins concentrate in specific regions of the arity map? In particular, do they overlap with the 7×13 hot spot the paper flags as anomalous in Q3?
- This extension is pending final confirmation with the professor.

**End of Week 2 milestone:** all proteome-scale figures reproduced; results table ready.

### Week 3 — Writeup and polish

**Days 15–17: LaTeX report**
- Standard structure: Introduction, Background (TDA primer, filtrations, persistence), Methods, Results, Discussion, References.
- All figures generated from the pipeline (no manual editing).
- Include the null-model comparison explicitly.

**Days 18–19: Final presentation slides**
- ~15 slides: motivation → AlphaFold and pLDDT → the filtration idea → the five statistics → null model → proteome results → fragmented subset → discussion.
- Reuse plots from the report.

**Day 20: Polish, sanity checks, buffer**
- Re-run the full pipeline from a clean checkout to verify reproducibility.
- Spot-check 5 random proteins by hand.
- Final pass on the report and slides.

**Day 21: Buffer / submission.**

---

## 6. Repository structure

```
project_q4/
├── README.md
├── requirements.txt
├── data/
│   ├── prototypes/                  # the 3 prototype PDBs
│   └── hsapiens/                    # full proteome (gitignored; ~11 GB)
├── src/
│   ├── parse.py                     # PDB/mmCIF → (sequence, pLDDT array)
│   ├── filtration.py                # Union-Find, PD, N_cc curve
│   ├── statistics.py                # f+_cp, p̄, H_p, N_cc^max, PLM
│   ├── null_model.py                # random pLDDT baselines, Conjecture 1
│   ├── plm.py                       # gudhi wrapper for 2nd-layer persistence
│   ├── process_proteome.py          # batch pipeline
│   └── plots.py                     # Figures 3, 7, 8
├── notebooks/
│   ├── 01_single_protein_walkthrough.ipynb
│   ├── 02_null_model.ipynb
│   ├── 03_prototypes_figure3.ipynb
│   └── 04_proteome_figures7_8.ipynb
├── results/
│   ├── proteome_stats.parquet
│   └── figures/
├── tests/
│   └── test_*.py
├── report/
│   ├── main.tex
│   ├── refs.bib
│   └── figures/                     # symlinked from results/figures
└── slides/
    ├── intro.pptx                   # Week 1 introductory deck
    └── final.pptx                   # Week 3 final presentation
```

---

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Proteome download is large (~11 GB) and slow | Pipeline supports a `--subset N` flag to run on a random sample first; the full run can happen overnight in Week 2. |
| pLDDT field convention mismatch | AlphaFold writes pLDDT in the B-factor column on a 0–100 scale; verify on the prototype proteins before batch-running. |
| Fragmented entries (F1, F2, …) inflate the protein count | Decide policy in Week 2 day 8 and document it; default = treat each fragment as a separate "protein" (matches the paper's protein count for H. Sapiens). |
| Morse-Smale simplification edge cases | Use gudhi's cubical complex on the N_cc array; well-tested. |
| Discrepancies with paper's exact numbers | Expected — paper's exact counts depend on AlphaFold-DB version. Report should explain version used and any deltas. |
| Time overrun | Day 14 extension is descoped from the report to a one-page appendix if behind schedule. |

---

## 8. What is NOT in scope

To keep the project focused:
- No reproduction of Q2 (ECOD domains) or Q3 (IDR/DisProt/AIUPred).
- From Q1, only the arity signature and arity-map binning are reused — for the cross-correlation extension in Week 2, Day 14. The full Q1 analysis (hierarchical clustering, optimal transport distances) is out of scope.
- No AlphaFold inference — we use AlphaFold-DB precomputed predictions.
- No comparison to AlphaFold 3 — paper uses AlphaFold 2 / AlphaFold-DB v4.

---

## 9. Decisions confirmed with the professor

- **Proteome scope:** H. Sapiens analysis is sufficient. No additional proteomes required.
- **Optional twist:** The Q1-arity cross-correlation is the planned extension (pending final confirmation with the professor). This is reflected in Section 5, Day 14.
- **Citation style:** `ieeetr` for the LaTeX report.
