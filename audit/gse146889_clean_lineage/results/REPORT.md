# GSE146889 clean-lineage confound-attribution rerun

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)

## Scope and provenance

This audit re-executed the GSE146889 confound-attribution analysis from the public GEO RPKM matrix under a protocol committed before execution. It did not import the historical transcript-reconstructed GSE146889 driver, did not load its saved null pickles, and computed all new statistics before reading historical result tables. It is an internal computational replication, not independent human review.

- Input SHA-256: `16e22f1285cb3960cc30151e67f8bf1ccc94a5d2e42051de9cc47ca1bbba3fa4`
- Input bytes: `53882628`
- Samples: 176 (91 tumor, 85 normal)
- Valid gene rows: 64253; invalid/missing rows removed: 1
- Protocol anchor commit: `5bedc6c0f1e55764e2b8683eed91a975ed25191f`
- Execution commit: `2a3649d7a768001d61430e48689335d59c2d785f`

## Result

**Historical qualitative verdict sustained.**

Deterministic real-data values reproduced within an absolute tolerance of 0.05: **False** (maximum absolute difference `0.0869281`).

The clean rerun supports the following bounded interpretation:

- Aggregate max-H1 excess survives linear class-mean residualization against both nulls: **True**.
- All four confound-projection quartiles exceed z=3 against both nulls: **True**.
- Tumor-only and normal-only observed max-H1 values are materially below the mixed value: **True**.
- Both single-class row-bootstrap intervals include zero: **True**.

This sustains a qualified claim: the aggregate topological excess is not removed by the exact linear tumor/normal mean shift, but the single-class signal is weaker and unstable under row resampling. It does not establish complete confound independence, causality, or external biological validity.

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

| space | observed_max_H1 | z_gaussian | z_permutation | cv_auc_tumor_normal | top_loop_n_participants | top_loop_n_tumor | top_loop_n_normal | fisher_p_tumor_enrichment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full HVG->PCA50 (class-mean intact) | 7.5643 | 11.3926 | 11.7631 | 0.9494 | 45 | 40 | 5 | 1.9405e-09 |
| Class-mean-residualized HVG->PCA50 | 8.4149 | 13.6345 | 13.9849 | 0.1191 | 22 | 10 | 12 | 0.8037 |

The loop-composition columns use support vertices of ripser's representative H1 cocycle. They are not a canonical cycle identity and are reported only as a descriptive sensitivity diagnostic.

## Bootstrap stability

| subset | n | point_delta | delta_ci95_lo | delta_ci95_hi | excludes_zero |
| --- | --- | --- | --- | --- | --- |
| mixed (tumor+normal) | 176 | 4.3228 | 2.5583 | 9.0135 | yes |
| tumor-only | 91 | 2.3500 | -0.6389 | 4.0307 | no |
| normal-only | 85 | 2.8394 | -1.5095 | 4.1822 | no |

The bootstrap holds the Gaussian null distribution fixed across row resamples, matching the historical audit. This is an approximation and is not a joint bootstrap of both observed and null processes.

## Differences from the historical artifact

`historical_comparison.csv` contains every directly comparable value. Monte Carlo and bootstrap differences are expected because this rerun regenerated all random draws independently rather than reusing historical artifacts. Deterministic differences are the primary reproduction test.

| section | key | metric | metric_type | historical | new | absolute_difference | relative_difference | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| within_class | Tumor-only | p_gaussian | monte_carlo | 0.0033 | 0.0066 | 0.0033 | 1.0000 | different |
| residualization | Class-mean-residualized HVG->PCA50 | cv_auc_tumor_normal | deterministic | 0.2060 | 0.1191 | 0.0869 | 0.4219 | different |
| residualization | Class-mean-residualized HVG->PCA50 | fisher_p_tumor_enrichment | diagnostic | 0.6496 | 0.8037 | 0.1541 | 0.2373 | different |
| residualization | Full HVG->PCA50 (class-mean intact) | cv_auc_tumor_normal | deterministic | 0.9254 | 0.9494 | 0.0240 | 0.0259 | different |
| bootstrap | tumor-only | z_ci95_lo | bootstrap | -1.3928 | -1.2178 | 0.1749 | 0.1256 | different |
Only the first 30 differently classified entries are shown here; the CSV is authoritative.

## Data assessment

The public GEO supplementary RPKM matrix is a legitimate source for this cohort and the structural sample counts reproduce. The analysis remains limited by heterogeneous tissue composition (colorectal, endometrial, and ovarian samples), binary tumor/normal residualization that cannot remove all tissue and nonlinear structure, and data-dependent HVG/PCA geometry. Those limitations narrow interpretation but do not invalidate the computational rerun.

## Files

- `primary_results.json`: machine-readable conclusions and complete deterministic results
- `within_class.csv`, `within_stratum.csv`, `residualization.csv`, `bootstrap.csv`, `confound_spectrum.csv`: report tables
- `null_and_bootstrap_draws.npz`: every regenerated Monte Carlo and bootstrap draw
- `historical_comparison.csv`: new-vs-historical comparison without forced agreement
- `gse146889_clean_lineage_audit.png`: generated only from the new CSV outputs
- `run_provenance.json`: input, environment, script, shard, and output hashes
