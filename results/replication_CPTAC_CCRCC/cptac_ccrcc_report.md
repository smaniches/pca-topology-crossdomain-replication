# CPTAC CCRCC proteomics: pre-registered confirmatory test of PCA-driven persistent-homology inflation

**Author:** Claude Science (TOPOLOGICA project) | **Date:** 2026-07-08
**Governing document:** `PREREGISTRATION.md`, locked 2026-07-08T13:50:55Z, SHA-256 `5e53930974...`
**Status:** Confirmatory (dataset #4 of the pre-registration's planned program, Section 8)

## 1. Data access

**No substitution was needed.** CPTAC proteomics data was accessed directly via the public
Proteomic Data Commons (PDC) GraphQL API (`pdc.cancer.gov/graphql`) once network access to that
domain was granted. The pre-registration's contingency plan (Section 8: substitute a GEO/PRIDE
dataset if CPTAC access is blocked) was not triggered.

**Dataset:** CPTAC Clear Cell Renal Cell Carcinoma (CCRCC) proteome study, `PDC000127`
(TMT10, JHU/PNNL processing), retrieved via the `quantDataMatrix` API endpoint
(`data_type="log2_ratio"`). This returns **9,591 proteins × 194 samples** (110 Primary Tumor,
84 Solid Tissue Normal — matched via the study's `biospecimenPerStudy` metadata; 5 cell-line
and 9 not-reported aliquots were excluded from the 208-sample raw pull). Missing values
(6.5% of entries) are PDC's own non-detected/masked-value encoding.

## 2. Deviation from Section 4 (disclosed, per Section 4's own instruction)

Section 4 specifies "median-normalize across samples, then log2 transform" for CPTAC proteomics,
written assuming **raw linear-scale intensities** as input. The PDC public API, however, only
exposes CCRCC protein abundance as **log2(ratio-to-common-reference)** — already log-transformed,
with negative values throughout (range: −16.14 to 11.72). Applying a second log2 transform to
signed log-ratio data is undefined/not meaningful.

**Deviation applied:** (1) impute missing → 0 (Section 3's imputation rule, unchanged), (2)
median-**center** each sample in log-ratio space (the log-space equivalent of "median-normalize"),
(3) skip the redundant second log2 step since the data already arrives log-transformed. All
downstream steps (HVG-selection(2000)/standardize/PCA(50)) are identical to Sections 3–4,
unchanged. Observed as a side effect: PDC's log2_ratio matrix arrives from CPTAC's own upstream
TMT pipeline already per-sample median-centered (every sample's median log-ratio was exactly
0.0 pre-normalization) — our step is therefore a confirming no-op on this dataset, not a new
transformation.

## 3. Intrinsic-dimension gate (Section 5)

| Space | d_int (MLE) | d_int (TwoNN) | d_ambient | ratio | Regime |
|---|---|---|---|---|---|
| Raw top-2000 HVG (standardized) | 18.61 | 17.38 | 2000 | 0.0093 | JL |
| Raw bottom-2000 variance (standardized) | 20.07 | 15.86 | 2000 | 0.0100 | JL |
| PCA(50) | 10.22 | 9.58 | 50 | 0.2043 | JL |

All three conditions fall inside the Johnson-Lindenstrauss regime (ratio < 0.3). **No metric
cascade deviation was required** — standard Euclidean-distance Vietoris-Rips persistent homology
is methodologically valid throughout, exactly as in the GSE81089 pilot.

## 4. Real-data persistent homology (H0/H1/H2)

| Space | H0 (components) | H1 (loops), max persistence | H2 (voids), max persistence |
|---|---|---|---|
| Raw top-2000 HVG | 193 | 78 bars, **3.833** | 13 bars, 0.631 |
| Raw bottom-2000 variance | 193 | 50 bars, 3.410 | 10 bars, 0.910 |
| PCA(50) | 193 | 76 bars, **4.683** | 6 bars, 0.377 |

Unlike the GSE81089 pilot (where bottom-variance genes showed **zero** H1 structure, an
underpowered-null result attributed to bulk RNA-seq's low sample count), CPTAC's bottom-2000
lowest-variance proteins **do** carry detectable H1 topology (50 bars, max persistence 3.41) —
closer to VR-003's original scRNA finding that low-variance features carry real signal, though
this dataset's primary statistic is still dominated by the top-variance and PCA spaces.

## 5. Null models (both pipeline-symmetric, per Section 3 hyperparameters)

- **Gaussian null:** i.i.d. N(0,1), dimension-matched (194×2000), 500 independent draws, each
  passed through the identical standardize→PCA(50) pipeline.
- **Pipeline-symmetric permutation null:** each protein's values independently permuted across
  the 194 samples (preserves each protein's marginal distribution, destroys cross-protein/
  cross-sample correlation structure), 2000 permutations, identical downstream pipeline.

All stochastic steps used seed 42 as specified. Both null computations ran to completion
(500/500 Gaussian draws in 27 min; 2000/2000 permutations in 98 min).

## 6. Primary statistic results: max H1 persistence, real vs. null

| Comparison | Observed | Null mean ± SD | Raw p | BH-adjusted p | z-score (standardized deviation from null) |
|---|---|---|---|---|---|
| Real raw-HVG vs. Gaussian null | 3.833 | 1.030 ± 0.121 | 0.0020 | 0.0023 | **23.18** |
| Real PCA50 vs. Gaussian null | 4.683 | 3.126 ± 0.354 | 0.0040 | 0.0040 | **4.40** |
| Real raw-HVG vs. pipeline-symmetric null | 3.833 | 0.983 ± 0.133 | 0.0005 | 0.0010 | **21.39** |
| Real PCA50 vs. pipeline-symmetric null | 4.683 | 2.992 ± 0.364 | 0.0005 | 0.0010 | **4.64** |

**z-score labeling note:** as specified in the pre-registration, these are standardized
deviations from a null distribution — (observed − null mean) / null SD computed against a
single real-data observation — **not Cohen's d**, which would require variance on both sides of
a paired two-sample comparison. This label was corrected in the GSE81089 pilot's own report
after an earlier mislabeling; the same correct labeling is used here from the start.

BH-FDR correction (α=0.05) was applied across the full pre-registered family of up to 20 tests
(Section 6): the 4 GSE81089 pilot tests plus these 4 CPTAC tests (8 of the ≤20-test family
populated so far). All 8 raw p-values and BH-adjusted p-values are reported in
`bh_correction_family_table.csv`; **all 8 remain significant after correction.**

## 7. Section 7 effect-size threshold: H0 verdict

The pre-registered bar requires z > 3.0 in **all four** required comparisons (raw × {Gaussian,
permutation}, PCA × {Gaussian, permutation}).

| Comparison | z-score | Pass (z > 3.0)? |
|---|---|---|
| Raw HVG vs. Gaussian null | 23.18 | ✅ |
| PCA50 vs. Gaussian null | 4.40 | ✅ |
| Raw HVG vs. pipeline-symmetric null | 21.39 | ✅ |
| PCA50 vs. pipeline-symmetric null | 4.64 | ✅ |

**H0 VERDICT: PASS.** CPTAC CCRCC replicates the original signal-vs-null effect — real
max-H1-persistence significantly exceeds both null models, in both raw and PCA-reduced feature
space, clearing the pre-registered z > 3.0 bar in all four required comparisons (weakest margin:
z=4.40, PCA vs. Gaussian — comparable to the pilot's own weakest margin of z=5.2, PCA vs.
pipeline-null).

## 8. Section 1 H1 hypothesis: does PCA inflation exceed a matched null?

| Quantity | Value |
|---|---|
| Real PCA-delta (PCA50 − raw) | **+0.850** |
| Gaussian-null PCA-delta (mean ± SD, paired per-draw) | **+2.096 ± 0.372** |
| Pipeline-null PCA-delta (mean ± SD, paired per-perm) | **+2.009 ± 0.386** |
| Gaussian null ≥ real? | **True** |
| Pipeline null ≥ real? | **True** |

Paired Wilcoxon signed-rank tests confirm PCA inflates max-H1-persistence within both null
models (Gaussian: p=1.26e-83; pipeline-symmetric: p≈0), with paired Cohen's d of 5.64 and 5.21
respectively — both very large effects, legitimately labeled as Cohen's d here because this is
a genuine paired two-condition (PCA vs. raw) comparison within each null's own draws.

**H1 VERDICT: SUPPORTED.** Both null models show a larger PCA-driven persistence increase
(+2.0 to +2.1) than the real data's own PCA-driven increase (+0.85). Per Section 1's decision
rule ("a dataset is judged FOR this hypothesis if null PCA-delta ≥ real PCA-delta, within each
null model separately"), CPTAC CCRCC supports the pre-registered H1: **PCA-driven persistence
increase is not, by itself, evidence of PCA revealing more real structure — a matched null
(pure noise or permuted data) shows an equal-or-larger increase from the same geometric
projection artifact.** This replicates the pilot's directional finding (GSE81089: real=+0.78,
Gaussian-null=+2.00, pipeline-null=+1.50) almost exactly in magnitude and direction on an
independent cancer type and an independent omics modality (proteomics vs. transcriptomics).

## 9. Confound check

**Classification confound (CV AUC):** 5-fold stratified cross-validated logistic regression on
PCA(50) coordinates predicts tumor/normal status with **AUC = 1.000 ± 0.000** — as expected,
CCRCC tumor vs. normal proteomes are near-perfectly separable, so overall sample structure in
this dataset is dominated by the tumor/normal axis.

**Cocycle trace (which samples does the top H1 loop touch?):**

| Space | Max H1 persistence | Loop vertices (n) | Tumor/Normal composition of loop | Fisher's exact p |
|---|---|---|---|---|
| Raw top-2000 HVG | 3.833 | 3 | 3 tumor, 0 normal | 0.260 |
| PCA(50) | 4.683 | 4 | 4 tumor, 0 normal | 0.135 |
| Bottom-2000 variance | 3.410 | 5 | 4 tumor, 1 normal | 0.391 |

None of the three top H1 loops show a statistically significant tumor/normal label enrichment
(all Fisher's exact p > 0.13). As in the GSE81089 pilot, **the most persistent H1 loop is not
simply encoding the tumor/normal label** — but with only 3–5 participating samples out of 194,
this check is underpowered to rule out a narrower tumor subtype or batch/technical origin. All
three top loops happen to skew toward all-tumor or majority-tumor composition, which is
suggestive but not statistically confirmed at this scale.

## 10. Honest conclusion

**H0 (does the original signal-vs-null effect replicate): PASS.** CPTAC CCRCC clears every
pre-registered effect-size threshold. Real proteomic data carries genuine topological (H1)
structure exceeding both a Gaussian-noise null and a pipeline-symmetric permutation null, in
both raw and PCA-reduced space — a second independent confirmation of "real omics data carries
real topology beyond a matched null," this time in a different omics modality (protein
abundance, not transcript abundance) and a different cancer type (renal, not lung).

**H1 (is PCA inflation a null-comparable artifact, not evidence of more real structure):
SUPPORTED.** The real data's PCA-driven persistence increase (+0.85) is smaller than the
artifactual inflation seen in both null models (+2.0–2.1) — closely replicating the pilot's
core methodological finding. This is now the second independent dataset (different modality,
different cancer type) showing that naively interpreting a PCA-space topological increase as
"PCA revealing more real structure" is unsupported without a matched null, because pure or
permuted noise shows an equal-or-larger increase from the same projection geometry.

**Confound check:** the top H1 loop is not a simple proxy for the tumor/normal label
(Fisher's exact p > 0.13 in all three feature spaces), though the check is underpowered
(3–5 participating samples) to determine what the loop does represent.

## 11. Limitations

- **Single confirmatory dataset per modality.** This is one CPTAC cohort (CCRCC); other CPTAC
  cancer types were not tested and generalization within proteomics is untested beyond this one
  cohort.
- **log-ratio, not linear-intensity, input.** The Section 4 normalization rule was written for
  raw linear intensities; the disclosed deviation (median-centering log-ratios, skipping the
  redundant log2) is the closest-possible adherence given what PDC's public API actually exposes,
  but is not byte-identical to the rule as literally written.
- **Missing-value imputation to 0.** 6.5% of entries were non-detected/masked and imputed to 0
  per Section 3's rule; this is a standard but lossy convention for proteomic missingness (which
  is often not missing-at-random — low-abundance proteins are systematically more likely to be
  undetected).
- **Confound check underpowered.** As in the pilot, only 3–5 samples participate in the top H1
  loop, too few for a well-powered Fisher's exact test; a larger cohort or a loop-membership
  bootstrap would strengthen this check.
- **BH-correction family is still accumulating.** Per Section 6, the full family spans up to
  20 tests (5 datasets × 2 nulls × 2 spaces); this analysis populates 8 of those 20. The
  correction reported here (family of 8) will be revised as further pre-registered datasets are
  analyzed, per Section 6's own accounting rule.

## 12. Artifacts

- `null_distributions.png` — null distributions (Gaussian, pipeline-symmetric) vs. real observed max-H1-persistence, raw and PCA50 spaces
- `pca_delta_and_persistence_diagram.png` — PCA-delta comparison (H1 hypothesis test) and real-data H1 persistence diagram (raw HVG space)
- `cptac_final_results_table.csv` — full H0 test statistics (observed, null mean/SD, z-score, p-values, BH-adjusted p, PASS/FAIL)
- `cptac_h1_pca_delta_table.csv` — H1 hypothesis test values (real vs. null PCA-deltas)
- `bh_correction_family_table.csv` — full 8-test BH-FDR correction family (pilot + CPTAC)
