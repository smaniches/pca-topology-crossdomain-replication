# CPTAC-CCRCC class-conditional matched-null sensitivity (2026-10-09)

**Status: completed prospective sensitivity experiment (NOT preregistered confirmation).**
This report does **not** alter the July 2026 preregistration, original confirmatory tests, or historical confound-audit results.

## Why this was necessary

The historical confound-audit residualization control compares class-mean-residualized observed features with null draws for which the same class-mean residualization was not re-applied. The reported residualized-space cross-validated AUC also processes held-out observations using their true labels. Neither the previous residualization-null z-score nor the globally residualized AUC permutation p-value should be interpreted as independently validated biological-confound removal.

This experiment corrects *one* asymmetry: both observed and null features are residualized by class means, PCA(50)-transformed, and scored for max H1 using the same algorithm and random seed.

## Frozen analysis contract

- Repository base: `70f7f6db91f59839eadc66bcbc1151e03e248edc`.
- Executed research commit: `e48f86a32994dcdde7694f01ea638c1fd60c123a`; the 499-draw job succeeded as GitHub Actions run [37946929016](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/37946929016).
- Input: committed PDC000127 CPTAC-CCRCC proteomics cache, 194 samples (110 primary tumor / 84 normal), top 2,000 variable proteins under the existing cohort's impute-zero / per-sample median-centering / standardization logic.
- Input matrix SHA-256: `d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3`.
- Labels SHA-256: `5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e`.
- Null: for each draw, **independently permute every feature within tumor and normal strata separately**, preserving each observed class-conditional univariate empirical distribution but removing within-class cross-protein correlations. Leave labels fixed. Recompute class-mean residualization, PCA(50) with seed 42, and ripser max-H1 for every draw, identical to the observed-statistic processing.
- Null RNG: NumPy `default_rng(20261009)`; 499 draws; one-sided plus-one Monte Carlo p-value `(1 + count(null >= observed))/(499 + 1)`.
- Pinned analysis environment: Python 3.11.17; numpy 2.4.6; pandas 3.0.3; scipy 1.17.1; scikit-learn 1.9.0; ripser 0.6.14.
- Six synthetic/invariance tests passed. The 24-draw pilot passed first, and the 499-draw job was launched only after that evidence was reviewed.

## Results

| Metric | Value |
|---|---:|
| Observed class-residualized PCA(50) max-H1 | 4.4099388123 |
| Historical Gudhi max-H1 on this observed point cloud | 4.4099384252 |
| Matched conditional null mean | 0.8230674233 |
| Matched conditional null SD (sample) | 0.1306607076 |
| Matched conditional null minimum / maximum | 0.4603900909 / 1.2946519852 |
| Exceedances among 499 draws | 0 |
| One-sided Monte Carlo plus-one p-value | 0.002 |
| Descriptive (observed - null mean)/null SD | 27.4518 |
| Difference, observed minus null mean | 3.5868713889 |
| Scientific computation time, excluding environment install | 18.18 seconds |

The observed numerical statistic agrees with the historical Gudhi result to within 4e-7, now using ripser on both observed and null. Note that the historical null z-score and this conditional-null z-score have different reference distributions and **must not be compared as interchangeable significance estimates**.

## Interpretation and limitations

The data exceed the null distribution of independently permuted features **conditional on the two observed class labels**. This is compatible with multivariate dependence surviving class-mean removal, including ordinary *covariance* structure. It does **not** demonstrate nontrivial biological topology, independence from tumor/normal confounding, or specificity to persistent homology beyond simpler correlation geometry.

The null assumes that samples are exchangeable *within* each label stratum for each protein. We have not independently verified patient-pairing, batch effects, dependence among aliquots, or covariates from raw PDC metadata. The null destroys within-class cross-protein covariance and possible linked-sample structure; those are important alternative explanations. The next discriminating control is a **within-class covariance-spectrum-preserving null** with explicit batch/patient-dependence treatment, not an indiscriminate rerun of prior experiments.

No old p-values have been retrospectively corrected by this sensitivity run. The existing globally class-label-residualized cross-validation AUC permutation inference needs separate methodological reconciliation.

## Evidence and repeatability

- Exact 499 unrounded null draw values and all parameters: [499_draws.json](499_draws.json).
- The raw run artifact ZIP SHA-256 is `30e68dac5e5c40e0a191dd91d2ca66430d5687b3f630b968b5a47ddf6ab08f1f` (Actions artifact ID `11624447259`).
- The JSON file inside the artifact SHA-256 is `0d2410692ef799cee66921da1086ab32565ef4dd9272b02ddfda27e2e9601f31`.
- Source: [experiment.py](../../audit/cptac_resid_null/experiment.py); [six tests](../../audit/cptac_resid_null/test_experiment.py); [execution gate](../../audit/CPTAC_RESIDUALIZATION_NULL_GATE_20261009.md).
- The original results remain unchanged. The GitHub workflow returns to manual dispatch for 499-draw repetitions, with PR-triggered 24-draw pilots only.
