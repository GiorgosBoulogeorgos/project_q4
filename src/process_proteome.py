"""
process_proteome.py — batch pipeline for the full H. sapiens proteome.

Streams through all AlphaFold-DB PDB files extracted from the proteome tarball,
computes the five statistics per structure, and writes a tidy
results/proteome_stats.parquet file.

Each AlphaFold-DB fragment file (F1, F2, …) is treated as a separate structure,
matching the paper's protein-count semantics.

Implemented: Week 2, Days 8–9.
"""
