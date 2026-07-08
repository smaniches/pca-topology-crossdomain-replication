# Pre-registered confirmatory test #2: does PCA inflate persistent-homology signal beyond a matched null? GSE146889 (colorectal/endometrial/ovarian mismatch-repair-deficiency cohort)

**Author:** Claude Science (TOPOLOGICA project) | **Date:** 2026-07-08 | **Dataset:** GEO accession [GSE146889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146889) | **Pre-registration:** locked document, SHA-256 `5e539309747188a2e77aa36bf1f9aecd3b50b8493803efa756ee4715773f0517`, timestamp 2026-07-08T13:50:55Z

## 1. Scope and status

This is a **pre-registered confirmatory test**, run under the locked protocol governing all new dataset analyses in this program (the GSE81089 NSCLC pilot is excluded from this pre-registration — it is the prior that motivated it). No hyperparameter, threshold, or hypothesis wording was altered after the document was locked. All decisions below follow Sections 1-9 of the pre-registration exactly.

**Primary hypotheses under test:**
- **H0 (does the original signal-vs-null effect replicate):** real max-H1-persistence significantly exceeds both null models, in both raw and PCA-reduced space, at the pre-registered effect-size threshold (z > 3.0 in all four required comparisons).
- **H1 (is PCA-driven persistence inflation itself evidence of real structure, or a null-comparable artifact):** null PCA-delta ≥ real PCA-delta, within each null model separately.

## 2. Dataset selection and preprocessing

**Dataset:** GSE146889 — bulk RNA-seq (gene-level RPKM) from a mismatch-repair-deficiency study spanning **176 samples (91 tumor, 85 paired-normal)**, drawn predominantly from colorectal tissue with smaller endometrial and ovarian components (source-tissue breakdown from the series metadata: 80 colorectal, 72 endometrial, 2 ovarian, 22 unlabeled). This is distinct in cancer type/tissue from the GSE81089 pilot (NSCLC lung), satisfying the pre-registration's requirement for tissue-type independence. 64,254 annotated genes were provided in the raw count/RPKM matrix; one trailing footer row with no gene identifier was dropped, leaving 64,253 genes.

**Point-cloud convention:** samples are points (176 points), genes are ambient coordinates — identical convention to the pilot and to VR-003 (scRNA-seq).

**Preprocessing (Section 3-4, bulk RNA-seq rule — identical to pilot):** log1p(RPKM) transform; genes with zero variance dropped (51,634 of 64,253 retained); top-2000 highly-variable genes (HVG) selected by variance, and a separate bottom-2000-variance gene set retained; features standardized (zero mean, unit variance) before PCA; PCA computed to 50 components (88.7% variance explained on real data). Random seed 42 throughout.

**Imputation rule:** no negative/sentinel values were present in this RPKM matrix (unlike GSE81089's 0.008% `-1.0` sentinels); the imputation rule (missing → 0) had no cells to act on beyond the one dropped footer row, disclosed above.

**A caveat on tissue homogeneity:** unlike GSE81089 (single NSCLC cohort), this dataset mixes tissue types (colorectal/endometrial/ovarian) unified by a shared mismatch-repair-deficiency study design rather than a single organ. This is disclosed as a deviation from a strictly single-tissue design; it does not affect the pipeline's mechanics but should inform how "cross-cancer-type replication" is read for this dataset specifically.

## 3. Intrinsic-dimension gate (Section 5, metric-space-diagnostics protocol)

| Space | d_int (MLE) | d_int (TwoNN) | d_ambient | ratio (MLE) | ratio (TwoNN) | Regime |
|---|---|---|---|---|---|---|
| Raw HVG (top-2000) | 13.99 | 17.01 | 2000 | 0.007 | 0.009 | JL |
| Bottom-variance genes (2000) | 32.06 | 2.43 | 2000 | 0.016 | 0.001 | JL |
| PCA(50) | 8.53 | 8.97 | 50 | 0.171 | 0.179 | JL |

All three conditions fall well inside the Johnson-Lindenstrauss regime (ratio < 0.3 by both estimators), so standard Euclidean-distance Vietoris-Rips persistent homology is methodologically valid throughout. **No deviation from Euclidean distance was required** — no spectral or Fermat-distance fallback was applied, matching the pilot's finding.

## 4. Two independent null models (Section 3, identical hyperparameters to pilot)

1. **Gaussian null** — i.i.d. N(0,1) noise, dimension-matched (176×2000), independently re-drawn **500 times**, each through the identical standardize→PCA(50) pipeline as the real data.
2. **Pipeline-symmetric permutation null** — each gene's expression values independently permuted across the 176 samples, then passed through the identical HVG-reselection/standardization/PCA(50) pipeline, repeated **2000 times**.

## 5. Results

### 5.1 Topological summary (single-draw observed values)

| Condition | H0 (components) | H1 (loops) | H2 (voids) | Max H1 persistence |
|---|---|---|---|---|
| Real, raw HVG (2000 genes) | 176 | 75 | 18 | 6.204 |
| Real, bottom-variance genes (2000 genes) | 176 | **0** | 0 | 0.000 |
| Real, PCA(50) | 176 | 79 | 24 | 7.564 |

As in the pilot, the bottom-2000-variance genes show **zero detectable H1 topology** — this again does not replicate VR-003's scRNA finding (bottom-844 genes carrying the strongest H1 signal). With n=176, this bulk cohort remains far smaller than a typical scRNA dataset; the result is reported as observed, consistent with the pilot's own power caveat.

### 5.2 Permutation test: is observed max-H1-persistence higher than either null?

Primary statistic: **max H1 persistence**, per pre-registration Section 2.

| Comparison | Observed | Null mean ± SD | Raw p-value | z-score (standardized deviation from null) | Validation method |
|---|---|---|---|---|---|
| Real raw HVG vs. Gaussian null | 6.204 | 1.002 ± 0.116 | **0.0020** | **44.9** | Monte Carlo Gaussian null, n=500 |
| Real raw HVG vs. pipeline-symmetric null | 6.204 | 0.980 ± 0.132 | **0.0005** | **39.6** | Permutation test (per-gene shuffle), n=2000 |
| Real PCA50 vs. Gaussian null (PCA50) | 7.564 | 3.184 ± 0.365 | **0.0020** | **12.0** | Monte Carlo Gaussian null + PCA, n=500 |
| Real PCA50 vs. pipeline-symmetric null (PCA50) | 7.564 | 3.063 ± 0.396 | **0.0005** | **11.4** | Permutation test + identical PCA pipeline, n=2000 |
| Real bottom-variance genes vs. Gaussian null | 0.000 | 1.002 ± 0.116 | 1.0000 | -8.6 | Monte Carlo Gaussian null, n=500 |
| Real bottom-variance genes vs. pipeline-symmetric null | 0.000 | 0.980 ± 0.132 | 1.0000 | -7.4 | Permutation test (per-gene shuffle), n=2000 |

**Note on z-score labeling (per pre-registration correction from the pilot):** the "z-score" column above is (observed − null mean)/null SD computed against a single real-data observation. This is a standardized deviation from a null distribution, **not** a two-sample Cohen's d.

**Both null models agree: the real dataset's most persistent H1 feature is far outside either null distribution**, in both raw and PCA-reduced space (z = 11.4–44.9), all raw p-values at or near the resolution floor of their respective null (p=0.0005–0.0020). Raw p-values are reported unadjusted per Section 8 of the pre-registration; BH correction across the full cross-dataset family is deferred to the program-level analysis (Section 6).

### 5.3 H0 verdict against the pre-registered effect-size threshold (Section 7)

The pre-registered threshold requires z > 3.0 in **both** raw and PCA space, against **both** null models.

| Comparison | z-score | Threshold (z > 3.0) met? |
|---|---|---|
| Raw HVG vs. Gaussian null | 44.9 | YES |
| Raw HVG vs. pipeline-symmetric null | 39.6 | YES |
| PCA50 vs. Gaussian null | 12.0 | YES |
| PCA50 vs. pipeline-symmetric null | 11.4 | YES |

**H0 verdict: PASS.** All four required comparisons clear the pre-registered z > 3.0 threshold, several by a wide margin (raw-space z = 39.6-44.9; PCA-space z = 11.4-12.0). Real bulk RNA-seq topology significantly exceeds both a Gaussian-noise null and a pipeline-symmetric permutation null.

### 5.4 Does PCA inflate persistence in noise too? (H1 judgment)

| Comparison | Mean Δ(PCA − raw) | Paired test | Cohen's d (paired) |
|---|---|---|---|
| Gaussian null: PCA50 vs raw (paired, same 500 draws) | +2.18 | Wilcoxon p=1.26e-83 | 5.60 |
| Pipeline-symmetric null: PCA50 vs raw (paired, same 2000 perms) | +2.08 | Wilcoxon p=0.00e+00 | 5.10 |
| **Real data: PCA50 vs raw** | **+1.36** | (single observation, no distribution) | — |

**H1 verdict: SUPPORTED, replicating the pilot's finding.** PCA-driven inflation of max-H1-persistence in *both* pure-noise null models (+2.18 Gaussian, +2.08 pipeline-symmetric) **exceeds** the PCA-driven increase observed in the real data itself (+1.36). Per the pre-registration's decision rule (Section 1), this dataset is judged **FOR** H1: PCA-driven persistence increase is not, by itself, evidence of real structure being revealed — a matched null shows the same or larger increase from pure or permuted noise. This is a second, independent confirmation of the pilot's OQ-006 finding (PCA fabricates phantom H1 inflation on noise), now replicated in a distinct tissue/cancer-type cohort.

### 5.5 Confound check (Section 9 protocol, identical to pilot)

A 5-fold cross-validated logistic regression on PCA(50) coordinates predicts tumor/normal status with **AUC = 0.925 ± 0.019** — high but not perfect, unlike the pilot's AUC=1.000, consistent with this cohort's more heterogeneous tissue composition and mismatch-repair-status stratification.

Tracing the most persistent H1 loop's participating vertices via `ripser`'s cocycle representatives: the loop touches **45 of 176 samples** (40 tumor, 5 normal). A Fisher's exact test finds this tumor/normal composition **significantly enriched for tumor samples** relative to the full cohort's tumor:normal ratio (odds ratio = 12.5, p = 2.87e-09). **Unlike the pilot** (where the loop touched only 3 samples with no significant label enrichment), here the loop's sample composition **is significantly associated with tumor status** — this most persistent H1 feature is more plausibly tracking (or correlated with) the tumor/normal distinction, or a variable correlated with it (e.g. MMR-deficiency subgroup), rather than being clearly independent of the primary label as in the pilot.

## 6. Honest conclusion

**H0 (does the original signal-vs-null effect replicate): PASS.** Real bulk RNA-seq data from an independent cohort, in a different cancer type/tissue from the pilot, again shows persistent-homology signal in H1 significantly exceeding both a Gaussian-noise null and a pipeline-symmetric permutation null (raw p ≤ 0.0005, z = 11.4–44.9), clearing the pre-registered z>3.0 threshold in all four required comparisons. This directly replicates the pilot's H0 finding in a second, independent, cross-tissue dataset.

**H1 (is PCA inflation itself evidence of real structure, or a null-comparable artifact): supported — PCA inflation is null-comparable, again.** The real data's raw→PCA change (+1.36) is smaller than the PCA-driven inflation measured in both null models (+2.08 to +2.18). This replicates the pilot's central methodological finding: **max-H1-persistence is not, by itself, a reliable indicator of whether PCA is revealing or fabricating structure** — the statistic must be interpreted relative to a matched null, and here, as in the pilot, PCA amplifies persistence in noise by an amount comparable to or exceeding its effect on the real data.

**Confound check: partially different from the pilot.** Unlike GSE81089 (where the top H1 loop touched only 3 samples with no significant label association), here the top loop touches 45 of 176 samples and is significantly enriched for tumor status (Fisher p=2.9e-09). This does not overturn the H0/H1 verdicts above (which concern the *overall* topological signal relative to null, not this specific loop's biological interpretation), but it means the biological interpretability of the single most-persistent loop differs across datasets — in this cohort it is more plausibly tracking tumor status (or a correlated variable such as MMR-deficiency subtype) rather than being independent of it.

**Bottom-variance-gene replication of VR-003: NOT REPLICATED (as in the pilot).** Zero detectable H1 structure in the bottom-2000-variance genes, consistent with the pilot's non-replication and likely reflecting the same bulk-vs-single-cell power gap noted there.

## 7. Limitations

- **Mixed tissue composition.** This cohort spans colorectal (majority), endometrial, and ovarian tissue unified by a shared MMR-deficiency study design, not a single-organ cohort like GSE81089. This is a disclosed deviation from tissue-homogeneity and should be weighed when assessing "cross-cancer-type replication" for this specific comparison.
- **Confound-check divergence from pilot.** The top H1 loop's association with tumor status here (unlike the pilot) means this dataset's topological signal is less cleanly separable from the tumor/normal label than in GSE81089; a fully independent topological signal (uncorrelated with the primary clinical label) is not established here.
- **Sample count / power.** n=176 remains far below typical scRNA-seq scale; bottom-variance-gene null results should be read as underpowered rather than a genuine absence of structure.
- **Single real-data observation.** As in the pilot, the real dataset provides one observed value per statistic (not a distribution), appropriate for a permutation test but precluding a confidence interval on the real-data statistic itself.
- **Multiple-testing correction deferred.** Per Section 6 of the pre-registration, BH correction is applied later across the full cross-dataset family (up to 20 tests); the p-values reported here are raw and should not be read as already FDR-controlled.

## 8. Artifacts

- `persistence_diagrams_GSE146889.png` — H1 persistence diagram (real, raw HVG vs PCA50) and null-distribution histograms with observed value marked
- `barcode_and_pca_delta_GSE146889.png` — top-30 H1 barcode (raw HVG) and PCA-delta comparison (real vs. both nulls)
- `final_results_table_GSE146889.csv` — full permutation test statistics table (p-values, z-scores, validation method per comparison)
- `diagrams_GSE146889.pkl` — raw persistence diagrams (H0/H1/H2) for all three real-data conditions
- `null_distributions_GSE146889.pkl` — full Gaussian (500) and pipeline-symmetric (2000) null draws, raw and PCA50 max-H1-persistence
