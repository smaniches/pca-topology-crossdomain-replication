# Derived-output reconciliation

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)  
**Scope:** GSE146889 clean-lineage confound-attribution audit  
**Status:** comparison-definition correction; no experiment rerun

The clean-lineage workflow correctly recomputed the raw-data, PCA, persistent-homology, null, residualization, cocycle-support, and bootstrap quantities. During review, two mismatched definitions were found only in the derived historical-comparison layer:

1. The historical report's tumor/normal CV AUC used logistic regression with `C=1.0`, but the first aggregate report compared it with the clean rerun's newly added `C=0.01` sensitivity.
2. The historical report used a two-sided Fisher exact test for cocycle-support composition, but the first aggregate report compared it with the clean rerun's one-sided enrichment test.

Both variants had already been computed from the same clean rerun. This reconciliation selects the matched historical definitions for reproduction claims and retains the alternatives as explicitly labelled sensitivity analyses.

No raw input, preprocessing choice, selected HVG set, PCA coordinates, persistence diagram, observed max-H1 value, Monte Carlo draw, bootstrap draw, or cocycle support was changed. The full experiment was not rerun.

After matched-definition reconciliation:

- all deterministic max-H1 values reproduce exactly;
- the historical `C=1.0` intact and residualized AUC values reproduce exactly;
- the historical two-sided Fisher p-values reproduce exactly;
- regenerated stochastic quantities preserve the same scientific decisions;
- the historical qualitative GSE146889 verdict remains sustained.

`reconcile_results.py` performs the correction from the committed machine-readable clean-run outputs and writes the corrected report, comparison table, residualization table, figure, provenance ledger, and `SHA256SUMS` registry.
