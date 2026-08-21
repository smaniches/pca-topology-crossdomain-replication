# Derived-output reconciliation

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)  
**Scope:** GSE146889 clean-lineage confound-attribution audit  
**Status:** comparison-definition correction; no experiment rerun

The clean-lineage workflow correctly recomputed the raw-data, PCA, persistent-homology, null, residualization, cocycle-support, and bootstrap quantities. During review, two mismatched definitions were found only in the derived historical-comparison layer:

1. The historical report's tumor/normal CV AUC used logistic regression with `C=1.0`, but the first aggregate report compared it with the clean rerun's newly added `C=0.01` sensitivity.
2. The historical report used a two-sided Fisher exact test for cocycle-support composition, but the first aggregate report compared it with the clean rerun's one-sided enrichment test.

Both variants had already been computed from the same clean rerun. The reconciled artifacts select the matched historical definitions for reproduction claims and retain the alternatives as explicitly labelled sensitivity analyses.

No raw input, preprocessing choice, selected HVG set, PCA coordinates, persistence diagram, observed max-H1 value, Monte Carlo draw, bootstrap draw, or cocycle support was changed. The full experiment was not rerun.

After matched-definition reconciliation:

- all deterministic max-H1 values reproduce exactly;
- the historical `C=1.0` intact and residualized AUC values reproduce exactly;
- the historical two-sided Fisher p-values reproduce exactly;
- regenerated stochastic quantities preserve the same scientific decisions;
- the historical qualitative GSE146889 verdict remains sustained.

`verify_reconciliation.py` is deliberately read-only. It independently reconstructs the matched historical AUC and Fisher definitions from the committed clean-run machine-readable object, checks deterministic max-H1 values against the historical tables, verifies the bootstrap and residualization decision boundaries, and asserts that the committed reconciled comparison table and machine-readable claims agree. It does not rewrite any result.

The exact workflow used for the successful full clean-lineage execution is archived as `EXECUTION_WORKFLOW.yml` outside `.github/workflows/`; it is provenance, not an active automation. `results/SHA256SUMS` is generated only after the final reconciled artifacts are fixed and is verified by ordinary read-only CI.
