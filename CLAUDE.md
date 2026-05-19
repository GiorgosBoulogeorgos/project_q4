# CLAUDE.md

Context file for Claude Code sessions on this project. Read this first.

## Project context

This is a course project for "Algorithms in Structural Bioinformatics" reproducing
Q4 of Cazals & Sarti (2025, bioRxiv 2024.11.16.623929 v4) at the H. Sapiens
proteome scale, plus a Q1 cross-correlation extension.

**Authoritative documents:**
- `docs/project_plan.md` — the full 3-week plan. Follow this. If something is
  ambiguous, ask before deviating.
- `docs/papers/2024_11_16_623929v4_full.pdf` — main reference paper (v4, Apr 2025).
- `docs/papers/2024_11_16_623929_full.pdf` — older version (Jan 2025), sometimes
  clearer on the null model and Example 1.

## What's locked in (do not re-litigate)

- Question: Q4 only (pLDDT fragmentation), with a Q1 arity-map extension
  (confirmed with the professor).
- Scope: H. Sapiens proteome from AlphaFold-DB v4. No other organisms.
- Stack: Python 3.11+, numpy, scipy, gemmi (PDB parsing), pandas, matplotlib,
  seaborn, gudhi (1D persistence for PLM only), tqdm.
- The path-graph filtration G_pLDDT and its persistence diagram are implemented
  *from scratch* using Union-Find with the Elder Rule. The second-layer persistence
  for PLM uses gudhi's CubicalComplex. This split is intentional.
- Report citation style: ieeetr.
- Languages: English in code, report, and slides.

## Key reproduction targets

Three prototype proteins (Figure 3 of the paper):

| Type       | UniProt ID    | Length  | Expected H_p |
|------------|---------------|---------|--------------|
| Ordered    | P15121        | 316 aa  | ≈ 0.04       |
| Disordered | A0A0G2L439    | 449 aa  | ≈ 0.27       |
| Mixed      | Q9VQS4        | 781 aa  | ≈ 0.28       |

Proteome-scale targets:
- Figure 7: Pearson r(f⁺_cp, H_p) ≈ 0.97 across H. Sapiens proteome.
- Figure 8: filter (#a.a. ≥ 200, H_p ≥ 0.25, PLM ≥ 3) at t_p = 0.025 yields ~86
  structures in H. Sapiens.

## Conventions

- pLDDT is stored in the B-factor column of AlphaFold PDB files, on a 0–100
  scale. Verify on the prototype proteins before any batch run.
- The filtration processes residues by *decreasing* pLDDT (equivalent to inserting
  by increasing −pLDDT, which matches the paper's u = −pLDDT convention).
- Treat each AlphaFold-DB fragment file (F1, F2, …) as a separate structure.
  This matches the paper's protein-count semantics.
- Elder Rule: when two components merge, the *younger* one dies (gets a death
  pLDDT recorded in the PD); the older one persists.

## Workflow preferences

- TDD where it matters: every algorithm module (filtration, statistics, PLM)
  ships with a unit test covering the null random model and the three prototypes.
- Each logical step from `docs/project_plan.md` = one git commit (or a small
  PR-sized chunk for larger steps).
- After Week 1, before going proteome-scale, the three prototypes must reproduce
  the paper's H_p values within ±0.03.
- When you hit a check-in point in a prompt, stop and show output. Do not
  silently roll into the next step.

## Sanity checks to run early

- Inverted sort: if the disordered prototype gives H_p ≈ 0.04 instead of ~0.27
  (or vice versa), the pLDDT sort order is wrong. Should be **decreasing**.
- Null model: for n=1000 random pLDDT ∈ [0, 100], f⁺_cp should land near 0.33
  (paper's Example 1). If it's near 0.5 or 0.67, the accretion test is wrong —
  accretions are residues inserted with two already-present neighbours that get
  immediately swallowed (birth pLDDT == death pLDDT).
- pLDDT range: AlphaFold-DB writes pLDDT on a 0–100 scale, but some downstream
  tools (and some papers) use 0–1. Stay on 0–100 throughout to match the paper.

## Out of scope

- No Q2 (ECOD domains) or Q3 (DisProt/AIUPred).
- From Q1, only the arity signature and 2D arity-map binning (for the extension).
  The full Q1 analysis (hierarchical clustering, optimal transport distances) is
  not part of this project.
- No AlphaFold inference — we use AlphaFold-DB precomputed predictions only.
- No AlphaFold 3 comparison.
