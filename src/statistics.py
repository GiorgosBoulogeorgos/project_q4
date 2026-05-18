"""
statistics.py — five pLDDT fragmentation statistics.

Computes from the persistence diagram and N_cc curve:
    f+_cp   fraction of critical points with positive persistence
    p̄       mean persistence
    H_p     normalised Shannon entropy of the persistence distribution
    N_cc^max  peak connected-component count over the filtration sweep
    PLM(t_ν)  number of persistent local maxima of N_cc (see plm.py)

Implemented: Week 1, Days 3–4.
"""
