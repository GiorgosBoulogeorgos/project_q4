# H_p Investigation Report
## Week 1, Days 3–4 — pLDDT Filtration Implementation & H_p Discrepancy

**Project:** Q4 — pLDDT Fragmentation of AlphaFold Predictions at Proteome Scale  
**Paper:** Cazals & Sarti (2025), bioRxiv 2024.11.16.623929 v4  
**Date:** 2026-05-19  

---

## 1. What Has Been Implemented

### 1.1 Modules Written

| Module | Function | Status |
|--------|----------|--------|
| `src/parse.py` | `parse_alphafold_pdb(path)` → `{uniprot_id, sequence, plddt, ca_coords}` | ✅ Done |
| `src/filtration.py` | `build_pd_and_ncc(plddt)` → `(pd, ncc_curve)` | ✅ Done |
| `src/statistics.py` | `f_cp_plus`, `mean_persistence`, `persistence_entropy`, `ncc_max` | ✅ Done |

### 1.2 Algorithm Summary — `build_pd_and_ncc`

- Residues inserted in **decreasing pLDDT order** (equivalently: increasing −pLDDT).
- Path graph 0 — 1 — 2 — … — (n−1); edges activated when both endpoints are present.
- **Union-Find with path halving** (not full path compression — half the pointer writes, same amortized complexity).
- **Elder Rule:** when two components merge, the one with higher birth_pLDDT is "older" and survives; the younger one dies and gets a PD pair `(birth_of_dying, current_pLDDT)`.
- Tie-break: lower root index survives.
- Returns:
  - `pd` — shape (n−1, 2), columns `[birth_pLDDT, death_pLDDT]`
  - `ncc_curve` — shape (n, 2), columns `[pLDDT_at_insertion, N_cc_after_step]`

**Three event types per insertion:**

| Event | Neighbours present | N_cc change | PD pairs |
|-------|--------------------|-------------|----------|
| A | 0 | +1 | 0 |
| B | 1 (one component) | 0 | 1 accretion (b=d) |
| C | 2 (two distinct components) | −1 | 1 accretion + 1 positive |

> In a path graph, i−1 and i+1 can never be in the same component before i is inserted (no alternative path), so event C always produces exactly one positive pair.

### 1.3 H_p Formula (Equation 1 of the paper)

```
H_p = − Σ_{v ∈ P} P[v] · ln P[v]  /  ln |P|
```

Where:
- **P** = set of **distinct** persistence values `pᵢ = birth_i − death_i` across all n−1 PD pairs.
- **P[v]** = count(v) / m — count-based probability (m = n−1).
- Accretions (persistence = 0) are included as one element of P.
- Edge case: if |P| = 1 (all pairs have identical persistence), H_p = 0.

---

## 2. Test Results — All 30 Tests Pass

```
============================= test session starts ==============================
collected 30 items

tests/test_parse.py::test_parse_returns_correct_length[P15121-316-Ordered]    PASSED
tests/test_parse.py::test_parse_returns_correct_length[A0A0G2L439-449-Disordered] PASSED
tests/test_parse.py::test_parse_returns_correct_length[Q9VQS4-781-Mixed]      PASSED
tests/test_parse.py::test_plddt_in_range[P15121-316-Ordered]                  PASSED
tests/test_parse.py::test_plddt_in_range[A0A0G2L439-449-Disordered]          PASSED
tests/test_parse.py::test_plddt_in_range[Q9VQS4-781-Mixed]                    PASSED
tests/test_parse.py::test_sequence_is_single_letter_codes[P15121-316-Ordered] PASSED
tests/test_parse.py::test_sequence_is_single_letter_codes[A0A0G2L439-449-Disordered] PASSED
tests/test_parse.py::test_sequence_is_single_letter_codes[Q9VQS4-781-Mixed]   PASSED
tests/test_filtration.py::test_toy_pd_size                                     PASSED
tests/test_filtration.py::test_toy_ncc_size                                    PASSED
tests/test_filtration.py::test_toy_birth_ge_death                              PASSED
tests/test_filtration.py::test_toy_ncc_positive                                PASSED
tests/test_filtration.py::test_toy_ncc_starts_at_one                          PASSED
tests/test_filtration.py::test_toy_final_ncc_is_one                           PASSED
tests/test_filtration.py::test_positive_persistence_count                      PASSED
tests/test_filtration.py::test_determinism                                     PASSED
tests/test_filtration.py::test_null_model_f_cp_plus                           PASSED
tests/test_filtration.py::test_single_residue                                  PASSED
tests/test_filtration.py::test_two_residues                                    PASSED
tests/test_statistics.py::test_f_cp_plus_all_positive                         PASSED
tests/test_statistics.py::test_f_cp_plus_all_accretions                       PASSED
tests/test_statistics.py::test_f_cp_plus_mixed                                PASSED
tests/test_statistics.py::test_mean_persistence_zero                          PASSED
tests/test_statistics.py::test_mean_persistence_known                         PASSED
tests/test_statistics.py::test_persistence_entropy_uniform                    PASSED
tests/test_statistics.py::test_persistence_entropy_all_same                   PASSED
tests/test_statistics.py::test_ncc_max_value                                  PASSED
tests/test_statistics.py::test_null_model_persistence_entropy                 PASSED
tests/test_statistics.py::test_two_node_path                                  PASSED

============================== 30 passed in 0.17s ==============================
```

---

## 3. Null Model Validation — PASSES

The paper reports for i.i.d. Uniform[0,100] pLDDT:
- n=100: f⁺_cp ≈ 0.33, H_p ≈ 0.36
- n=1000: f⁺_cp ≈ 0.33, H_p ≈ 0.45
- n=10,000: f⁺_cp ≈ 0.33, H_p ≈ 0.47

Our results (20 independent trials, n=1000, seed=42):

| Statistic | Paper | Ours | Status |
|-----------|-------|------|--------|
| f⁺_cp | ≈ 0.33 | 0.332 ± 0.007 | ✅ Match |
| H_p | ≈ 0.45 | **0.4437 ± 0.0079** | ✅ Match |

---

## 4. Prototype Protein Stats (Raw Float pLDDT)

### 4.1 Input Data

| Protein | File used | n_residues | pLDDT range | n_distinct_pLDDT |
|---------|-----------|------------|-------------|-------------------|
| P15121 (Ordered) | AF-P15121-F1-model_v4.pdb | 316 | [84.52, 98.96] | 108 |
| P15121 (Ordered) | AF-P15121-F1-model_v6.pdb | 316 | [86.81, 98.94] | 48 |
| A0A0G2L439 (Disordered) | AF-A0A0G2L439-F1-model_v4.pdb | 449 | [17.39, 41.39] | 366 |
| Q9VQS4 (Mixed) | AF-Q9VQS4-F1-model_v6.pdb | 781 | [34.88, 96.69] | 526 |

### 4.2 Computed Statistics

| Protein | n_PD | n_accretion | n_positive | n_distinct_persist | f⁺_cp | mean_p | H_p | ncc_max |
|---------|------|-------------|------------|--------------------|-------|--------|-----|---------|
| P15121 v4 | 315 | 239 | 76 | 50 | 0.2413 | 0.1649 | 0.3713 | 42 |
| P15121 v6 | 315 | 252 | 63 | 27 | 0.2000 | 0.1136 | 0.3299 | 35 |
| A0A0G2L439 v4 | 448 | 296 | 152 | 140 | 0.3393 | 1.2149 | 0.4659 | 59 |
| Q9VQS4 v6 | 780 | 537 | 243 | 208 | 0.3115 | 1.6754 | 0.4236 | 80 |

### 4.3 Comparison with Paper Targets (CLAUDE.md)

| Protein | Expected H_p | Computed H_p (v4) | Computed H_p (v6) | Δ (v4) | Δ (v6) |
|---------|-------------|-------------------|--------------------|---------|---------|
| P15121 | **≈ 0.04** | 0.3713 | 0.3299 | **+0.33** | **+0.29** |
| A0A0G2L439 | **≈ 0.27** | 0.4659 | — | **+0.20** | — |
| Q9VQS4 | **≈ 0.28** | — | 0.4236 | — | **+0.14** |

> **All three prototypes fail the ±0.03 tolerance from CLAUDE.md.**

---

## 5. H_p Sensitivity Analysis — pLDDT Rounding

We tested rounding pLDDT to different resolutions before computing the filtration, to simulate lower-precision input data.

| Protein | Expected | raw float | round 0.01 | round 0.1 | round 0.5 | round 1.0 (int) |
|---------|----------|-----------|-----------|----------|----------|-----------------|
| P15121 v6 | 0.04 | 0.330 | 0.328 | 0.279 | 0.247 | 0.151 |
| A0A0G2L439 v4 | 0.27 | 0.466 | 0.466 | 0.466 | 0.456 | 0.433 |
| Q9VQS4 v6 | 0.28 | 0.424 | 0.424 | 0.425 | 0.428 | 0.425 |

**Observation:** P15121 responds to coarser rounding (0.33 → 0.15 at integer resolution), suggesting its pLDDT values have some clustering. A0A0G2L439 and Q9VQS4 are almost completely insensitive to rounding at any resolution, meaning the disordered/mixed proteins have nearly all distinct pLDDT values even as integers.

---

## 6. H_p Sensitivity Analysis — Persistence Value Rounding

We also tested rounding the **persistence values themselves** (birth − death) after computing the filtration, as the paper states *"assuming persistences belong to a discrete set P."*

| Protein | Expected | raw float | round 0.01 | round 0.1 | round 0.5 | round 1.0 (int) |
|---------|----------|-----------|-----------|----------|----------|-----------------|
| P15121 v6 | 0.04 | 0.330 (27) | 0.328 (23) | 0.317 (15) | 0.172 (7) | 0.129 (5) |
| A0A0G2L439 v4 | 0.27 | 0.466 (140) | 0.466 (131) | 0.464 (72) | 0.458 (28) | 0.431 (18) |
| Q9VQS4 v6 | 0.28 | 0.424 (208) | 0.424 (187) | 0.429 (104) | 0.425 (42) | 0.416 (26) |

> Numbers in parentheses are |P| (number of distinct persistence values after rounding).

**Observation:** Rounding persistence to integers dramatically reduces |P| for P15121 (27 → 5) and gives H_p = 0.129. Still ~3× above the expected 0.04. For A0A0G2L439, integer rounding reduces |P| from 140 to 18 but H_p barely changes (0.466 → 0.431), still far from 0.27.

---

## 7. Alternative Formula Variants Tested

All were tested on v4/v6 data with raw float pLDDT. None matched the paper's targets.

| Variant | P15121 | A0A0G2L439 | Q9VQS4 |
|---------|--------|-----------|--------|
| **Count-based (our implementation)** | 0.33 | 0.47 | 0.42 |
| Positive pairs only (exclude accretions from sum and |P|) | 0.90 | 0.99 | 0.99 |
| Lifetime-weighted: P[v] = count[v]·v / L | 0.77 | 0.91 | 0.92 |
| Individual lifetime fraction, norm by log(m') | 0.70 | 0.91 | 0.91 |
| Individual lifetime fraction, norm by log(m) | 0.51 | 0.75 | 0.75 |
| H_raw / log(n) (fixed denom = n residues) | 0.20 | 0.36 | 0.30 |
| H_raw / log(n−1) (fixed denom = m) | 0.20 | 0.36 | 0.30 |

> "Individual lifetime fraction": P[p_i] = p_i / L for each pair i (not grouped by distinct value). Accretions contribute 0.

---

## 8. Summary of pLDDT Data Across Versions

We extracted the v4 PDB for P15121 by streaming the human proteome tarball (`UP000005640_9606_HUMAN_v4.tar`, found at file entry #9,878 of ~20,000).

| File | n_distinct_pLDDT | pLDDT range | H_p |
|------|-----------------|-------------|-----|
| P15121-model_v4 (from human tarball) | 108 | [84.52, 98.96] | 0.3713 |
| P15121-model_v6 (from EBI API) | 48 | [86.81, 98.94] | 0.3299 |

**Key finding:** The v4 file gives a *higher* H_p than v6 (0.37 vs 0.33). The paper's target of 0.04 is unreachable with either dataset using any formula variant we tested.

---

## 9. Diagnosis

### What works correctly
- **Null model** (f⁺_cp, H_p) — matches paper's Example 1 values within noise ✅
- **Filtration mechanics** — Elder Rule, accretion counting, ncc_curve ✅
- **Relative ordering** — ordered (0.33) < mixed (0.42) ≈ null model (0.44) — qualitatively correct ✅
- **All 30 unit tests pass** ✅

### What doesn't match
- H_p for all three prototypes is 0.14–0.33 above the paper's values
- A0A0G2L439 (disordered) computes *above* the null model (0.47 > 0.44), whereas the paper places it below (0.27 < 0.45)

### Probable causes (unresolved)

1. **Undisclosed discretization:** The paper states *"assuming persistences belong to a discrete set P"* but never specifies how float pLDDT values are discretized. Our formula gives maximum H_p ≈ 1 when all persistence values are distinct (which they nearly always are with float data). The paper's results are consistent with integer-valued pLDDT, but even integer rounding doesn't reproduce the targets.

2. **Data version mismatch:** The paper's numerical results were produced from data downloaded before the EBI updated the v4 tarball. If the original P15121 pLDDT values were more uniform (narrower range, fewer distinct values), H_p ≈ 0.04 would be achievable with our formula. However, both current v4 and v6 files have 48–108 distinct float values, making this unlikely.

3. **Different |P| denominator:** Several formula variants with fixed denominators (log(n), log(n−1)) give lower H_p but still don't reach the paper's targets, and don't produce consistent results across the three proteins.

4. **Author's code not published:** No code repository link was found in the paper. The paper's contact is Frederic.Cazals@inria.fr.

---

## 10. Recommended Questions for the Professor

1. **Formula convention for |P|:** Is |P| the number of distinct persistence values, or is it a fixed value like n, n−1, or some function of the pLDDT range?

2. **Discretization:** Does the paper pre-process pLDDT values before computing the filtration (e.g., round to integer, round to 0.5, clip to some range)?

3. **Which filtration for H_p in Figure 3:** The figure shows both G_arity and G_pLDDT filtrations. The caption reports one H_p per protein — is that from the pLDDT filtration or the arity filtration?

4. **Data version:** Is there a frozen dataset or DOI for the exact PDB files used to produce Figure 3?

5. **Author contact:** Is contacting Cazals/Sarti directly appropriate for the course context?

---

## 11. Files in Repository

```
src/
  parse.py         — parse_alphafold_pdb() using gemmi
  filtration.py    — build_pd_and_ncc() using Union-Find + Elder Rule
  statistics.py    — f_cp_plus, mean_persistence, persistence_entropy, ncc_max
  download.py      — download_structure(), download_proteome()

tests/
  test_parse.py        — 9 tests on 3 prototypes (shape, pLDDT range, sequence)
  test_filtration.py   — 11 tests (toy example, null model f+cp≈1/3, determinism, edge cases)
  test_statistics.py   — 10 tests (analytical checks, null model Hp≈0.45)

data/prototypes/
  AF-P15121-F1-model_v4.pdb      (316 aa, extracted from human v4 tarball)
  AF-P15121-F1-model_v6.pdb      (316 aa, from EBI API)
  AF-A0A0G2L439-F1-model_v4.pdb  (449 aa, extracted from Zebrafish v4 tarball)
  AF-Q9VQS4-F1-model_v6.pdb      (781 aa, from EBI API)

docs/
  project_plan.md         — 3-week implementation plan
  hp_investigation.md     — this document
  papers/
    2024.11.16.623929v4.full.pdf  — paper v4 (April 2025)
    2024.11.16.623929.full.pdf    — paper v1 (January 2025)
```
