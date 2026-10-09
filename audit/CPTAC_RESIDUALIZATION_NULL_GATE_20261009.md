# CPTAC residualization-null research gate (2026-10-09)

Status: PROPOSED; no scientific result promoted.

Baseline: `70f7f6db91f59839eadc66bcbc1151e03e248edc`.
This is a prospective follow-up, not an amendment to the 2026-07-08 preregistration.

## Why this experiment
The existing `residualization_control` computes class-label-conditioned centering of the observed matrix, but its Gaussian and permutation references do not receive the same transformation. The subsequent AUC uses label-derived transformed features in cross-validation. Consequently, the AUC and post-residualization null comparisons cannot be treated as independently validated confound removal.

## Before any confirmatory compute
1. Freeze the exact estimand: max-H1 persistence of the class-mean-residualized CPTAC point cloud under a clearly specified metric and PCA fit.
2. Specify a conditional null model with the same label assignment and class-wise centering applied to every surrogate *before* fitting PCA and calculating max-H1. State explicitly what that null preserves and destroys; it is not a generic test of all confounding.
3. Verify the exchangeability assumptions for the proposed permutation. The class-mean residualization must be recomputed for each permuted matrix. A featurewise permutation may not preserve cross-feature covariance or within-patient dependence. Do not promote a p-value if exchangeability is unsupported.
4. Add deterministic unit tests: same-size outputs; zero class-mean difference; real and null transformation parity; fixed-seed reproducibility; positive/negative synthetic controls; PCA and H1 metric invariance checks.
5. Separate any label-predictability diagnostic from the topology test. Do not use globally label-conditioned residuals for out-of-fold AUC, and do not permute labels while holding label-derived residuals fixed.
6. Confirm the CPTAC committed input hashes from `checksums.sha256` and the current branch head.

## Compute gates
- Pilot: 20-32 draws on a standard public-repository hosted runner, with timeout <=30 minutes and no paid larger runner/GPU.
- Record measured wall time, memory, exceptions, surrogate-statistic distribution, and the exact software environment.
- Only if pilot and null-model validity checks pass, run 499 independent confirmatory draws using frozen code and seed, with a timeout <=120 minutes.
- A 499-draw Monte Carlo test has minimum attainable one-sided p-value 1/500 = 0.002; use the plus-one correction and report Monte Carlo uncertainty.
- Do not apply the historical z>3 rule as if this were preregistered confirmation; report it as a descriptive sensitivity alongside empirical tail probability and effect size.
- Stop on unexpected resource use, failed controls, unsupported null exchangeability, or missing source integrity.

## Outcomes
- If the corrected null retains excess observed max-H1, report evidence conditional on that null only.
- If not, weaken the residualization/confound-independence claim while preserving all historical results.
- Preserve the negative or inconclusive result; do not overwrite preregistered tables.
- Reconcile the manuscript's residualized-AUC permutation significance language with the GSE146889 clean-lineage reconciliation, which already withdrew one invalid p-value.

No claim of successful validation or compute execution is made by this document.
