"""
plm.py — gudhi wrapper for second-layer persistence (PLM statistic).

Feeds the 1D N_cc array into gudhi's CubicalComplex, extracts 0-dimensional
persistence pairs, and counts those with persistence ≥ t_ν = n · t_p.

This is the only module that uses gudhi; all other TDA is implemented from
scratch (see filtration.py).

Implemented: Week 1, Days 6–7.
"""
