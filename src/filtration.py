"""
filtration.py — path-graph filtration G_pLDDT and persistence diagram.

Builds the filtration by inserting residues in order of *decreasing* pLDDT
(equivalently: increasing −pLDDT, matching the paper's u = −pLDDT convention).
Uses Union-Find with the Elder Rule to record (birth, death) pairs and the
N_cc(pLDDT) curve.

Implemented: Week 1, Days 3–4.
"""
