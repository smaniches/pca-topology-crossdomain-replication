# GSE146889 clean-lineage confound-attribution rerun

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)

## Result

**Historical qualitative verdict sustained. Deterministic real-data quantities reproduced exactly under matched definitions.**

This audit re-executed the GSE146889 confound-attribution analysis from the public GEO RPKM matrix under a protocol committed before execution. It did not import the transcript-reconstructed driver, did not reuse its null pickles, and regenerated all Monte Carlo and bootstrap draws. It is an internal computational replication, not independent human validation.

- Input SHA-256: `16e22f1285cb3960cc30151e67f8bf1ccc94a5d2e42051de9cc47ca1bbba3fa4`
- Input bytes: `53882628`
- Samples: 176 (91 tumor, 85 normal)
- Valid gene rows: 64253; invalid/missing rows removed: 1

### Reconciliation of derived comparison semantics

The first aggregate report compared the historical AUC values (LogisticRegression `C=1.0`) against a newly added `C=0.01` sensitivity and compared historical two-sided Fisher exact p-values against one-sided enrichment p-values. Both variants were computed from the same clean rerun, but they are not like-for-like comparisons. This report uses the historical definitions for reproduction claims and retains the alternatives as explicitly labelled sensitivities. No raw-data, PH, null, or bootstrap result was changed or rerun.

- Deterministic maximum absolute difference under matched definitions: `0`.
- Aggregate residualized signal exceeds z=3 against both nulls: **True**.
- All four confound-projection quartiles exceed z=3 against both nulls: **True**.
- Both single-class bootstrap intervals include zero: **True**.

The bounded interpretation is unchanged: aggregate topological excess survives removal of the exact linear tumor/normal mean shift, while single-class signal is weaker and unstable under row resampling. This does not establish complete confound independence, causality, or external biological validity.

## Within-class decomposition

| condition | n | observed_max_H1 | z_gaussian | p_gaussian | z_permutation | p_permutation |
| --- | --- | --- | --- | --- | --- | --- |
| Mixed (tumor+normal) | 176 | 7.5643 | 11.3926 | 0.0020 | 11.7631 | 4.9975e-04 |
| Tumor-only | 91 | 5.8801 | 4.4792 | 0.0066 | 4.0871 | 0.0033 |
| Normal-only | 85 | 6.4021 | 5.4887 | 0.0033 | 5.2758 | 0.0033 |

## Within-stratum control

| quartile | n | n_tumor | n_normal | observed_max_H1 | z_gauss | z_perm |
| --- | --- | --- | --- | --- | --- | --- |
| Q1 | 44 | 1 | 43 | 1.9853 | 7.0557 | 6.6166 |
| Q2 | 44 | 9 | 35 | 2.7217 | 11.5337 | 9.9536 |
| Q3 | 44 | 37 | 7 | 2.3779 | 9.3128 | 8.2868 |
| Q4 | 44 | 44 | 0 | 4.4906 | 22.7722 | 20.2462 |

## Residualization

| space | observed_max_H1 | z_gaussian | z_permutation | cv_auc_tumor_normal | cv_auc_C_0_01_sensitivity | top_loop_n_participants | top_loop_n_tumor | top_loop_n_normal | fisher_p_tumor_enrichment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full HVG->PCA50 (class-mean intact) | 7.5643 | 11.3926 | 11.7631 | 0.9254 | 0.9494 | 45 | 40 | 5 | 2.8704e-09 |
| Class-mean-residualized HVG->PCA50 | 8.4149 | 13.6345 | 13.9849 | 0.2060 | 0.1191 | 22 | 10 | 12 | 0.6496 |

`cv_auc_tumor_normal` is the matched historical `C=1.0` result. `cv_auc_C_0_01_sensitivity` is retained as a regularization sensitivity. `fisher_p_tumor_enrichment` is the matched historical two-sided Fisher exact test; the one-sided value remains in `residualization.csv` as a labelled sensitivity. Cocycle support is a representative returned by ripser, not a canonical biological cycle.

## Bootstrap stability

| subset | n | point_delta | delta_ci95_lo | delta_ci95_hi | excludes_zero |
| --- | --- | --- | --- | --- | --- |
| mixed (tumor+normal) | 176 | 4.3228 | 2.5583 | 9.0135 | yes |
| tumor-only | 91 | 2.3500 | -0.6389 | 4.0307 | no |
| normal-only | 85 | 2.8394 | -1.5095 | 4.1822 | no |

The bootstrap holds the Gaussian null fixed across row resamples, matching the historical audit. This is not a joint bootstrap of observed and null processes.

## Differences from the historical artifact

Deterministic matched-definition quantities are exact. Remaining differences are regenerated Monte Carlo or bootstrap variation. `historical_comparison.csv` is authoritative.

| section | key | metric | metric_type | historical | new | absolute_difference | relative_difference | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| within_class | Tumor-only | p_gaussian | monte_carlo | 0.0033 | 0.0066 | 0.0033 | 1.0000 | different |
| bootstrap | tumor-only | z_ci95_lo | bootstrap | -1.3928 | -1.2178 | 0.1749 | 0.1256 | different |

## Data assessment

The public GEO supplementary RPKM matrix is the intended source and the structural cohort definition reproduced. The cohort mixes colorectal, endometrial, and ovarian tissue; binary tumor/normal residualization cannot remove all tissue, study-design, nonlinear, or latent biological structure. Those limitations constrain interpretation but do not invalidate the computational reproduction.

## Files

- `primary_results.json`: machine-readable results and corrected reproduction claims
- `within_class.csv`, `within_stratum.csv`, `residualization.csv`, `bootstrap.csv`, `confound_spectrum.csv`: report tables
- `null_and_bootstrap_draws.npz`: all regenerated stochastic draws
- `historical_comparison.csv`: matched-definition comparison
- `gse146889_clean_lineage_audit_reconciled.svg`: generated from corrected CSV outputs
- `gse146889_clean_lineage_audit.png`: preserved original workflow render before comparison-definition reconciliation
- `run_provenance.json`: immutable provenance record from the original clean-lineage workflow run
- `../RECONCILIATION.md`: derived comparison-definition correction, scope, and integrity rationale
- `SHA256SUMS`: audit-directory integrity registry
