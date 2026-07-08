# TCGA-LUAD Confirmatory Analysis Report
## Pre-registered test: Does PCA inflate persistent-homology signal beyond a matched null?

**Pre-registration:** locked 2026-07-08T13:50:55Z (SHA-256: 5e539309747188a2e77aa36bf1f9aecd3b50b8493803efa756ee4715773f0517)
**Dataset:** TCGA-LUAD (lung adenocarcinoma), fetched live via the GDC API (api.gdc.cancer.gov), Data Release 45.0.
**Status:** CONFIRMATORY (not exploratory) — all hyperparameters fixed before data was fetched.

---

## 1. Data provenance

- **RNA-seq**: STAR-Counts FPKM (`fpkm_unstranded`), GENCODE v36 gene model. 58 TCGA-LUAD cases with
  BOTH Primary Tumor and Solid Tissue Normal samples present (116 samples total, 60,660 genes).
  Files fetched individually via `GET /data/{file_id}` after a `POST /files` query for
  `data_category=Transcriptome Profiling`, `data_type=Gene Expression Quantification`,
  `analysis.workflow_type=STAR - Counts`.
- **DNA methylation**: SeSAMe Level-3 beta-values, Illumina Human Methylation 450K platform.
  18 of the 58 RNA-seq-paired cases also have paired tumor/normal methylation (36 samples,
  486,427 probes) — this 18-case subset is the SAME samples as a subset of the RNA-seq cohort,
  satisfying the task's "same samples" requirement for a genuine second omics layer.
- No synthetic or simulated data was used anywhere in this analysis. All matrices were built
  directly from GDC-downloaded files.

## 2. Post-hoc deviations (declared per Section 4, NOT silently absorbed)

**D1 — RNA-seq bottom-variance selection:** 5,556/60,660 genes have exactly zero variance across
all 116 samples (never detected), exceeding the 2,000-gene bottom-variance selection size. The
literal "bottom 2000 by variance" is therefore a degenerate all-identical (zero pairwise distance)
point cloud with trivial topology by construction. **Deviation applied:** bottom-2000 selected
among genes with variance > 0. Both the literal (degenerate, all-zero H0/H1/H2) and
deviation-applied results are reported.

**D2 — Methylation missing-value handling:** beta-values are bounded in [0,1]; the pre-registered
"missing → 0" rule is incompatible with `logit(beta)`, which diverges at beta=0. 64,886/486,427
probes are NaN in all 36 samples (dropped prior to HVG selection — no information to select on).
Remaining missing values are imputed as 0 (literal rule), then all values are clipped to
[1e-6, 1-1e-6] before the logit transform, applied uniformly (not only to imputed entries) for a
well-defined, symmetric transform.

**D3 — Methylation PCA components:** PCA(50) is infeasible with only 36 samples in the matched
methylation cohort (sklearn full-SVD cap = min(n_samples, n_features) = 36).
**Deviation applied:** PCA components reduced to `min(50, n_samples-1) = 35` for methylation only;
RNA-seq (116 samples) uses the pre-registered PCA(50) unmodified. A direct consequence: PCA(35) on
36 samples is a full-rank (100% variance retained), exactly isometric rotation of the standardized
raw space, so Euclidean-PCA max-H1 is numerically **identical** to the raw-space value — it is not
an independent test of PCA's effect and is reported for transparency only, not as a passing
comparison in the primary family below.

**D4 — Section 5 intrinsic-dimension gate mandated the metric cascade for methylation PCA:**
d_int/d_amb = 0.320 for the methylation PCA(35) space (MLE estimator), crossing the 0.3
"transitional" threshold. Per Section 5, Euclidean VR persistence is disallowed for this space and
the spectral (effective-resistance) distance was substituted. This is the metric used for all
methylation PCA-space primary-family comparisons below; the Euclidean-PCA numbers (identical to
raw per D3) are shown alongside only as a transparency check.

## 3. Intrinsic-dimension gate results (Section 5)

| Space | d_int (MLE) | d_int (TwoNN) | d_amb | ratio | regime | metric used |
|---|---|---|---|---|---|---|
| RNA-seq raw HVG (2000) | 14.5 | 16.0 | 2000 | 0.007 | JL | Euclidean |
| RNA-seq bottom-var (2000, deviation) | 22.8 | 71.4 | 2000 | 0.011 | JL | Euclidean |
| RNA-seq PCA(50) | 10.3 | 10.2 | 50 | 0.205 | JL | Euclidean |
| Methylation raw HVP (2000) | 11.2 | 16.1 | 2000 | 0.006 | JL | Euclidean |
| Methylation bottom-var (2000, deviation) | 17.5 | 15.2 | 2000 | 0.009 | JL | Euclidean |
| **Methylation PCA(35)** | **11.2** | **16.1** | **35** | **0.320** | **transitional** | **Spectral (cascade, D4)** |

## 4. Primary statistic: max H1 persistence, real vs. null (Section 7 test family)

BH-FDR correction (alpha=0.05) applied across this dataset's 8-test family (raw/PCA space x
Gaussian/permutation null x 2 omics layers), consistent with the up-to-20-test global family
defined in Section 6.

| test | real max-H1 | null mean | null std | z | p (raw) | p (BH) | reject @0.05 |
|---|---|---|---|---|---|---|---|
| RNAseq_raw_gaussian | 8.327 | 0.963 | 0.126 | 58.67 | 0.00200 | 0.00266 | True |
| RNAseq_raw_perm | 8.327 | 0.879 | 0.157 | 47.53 | 0.00050 | 0.00133 | True |
| RNAseq_pca50_gaussian | 7.433 | 3.470 | 0.431 | 9.19 | 0.00200 | 0.00266 | True |
| RNAseq_pca50_perm | 7.433 | 3.225 | 0.504 | 8.34 | 0.00050 | 0.00133 | True |
| Methylation_raw_gaussian | 4.092 | 0.793 | 0.171 | 19.34 | 0.00200 | 0.00266 | True |
| Methylation_raw_perm | 4.092 | 0.903 | 0.119 | 26.69 | 0.00050 | 0.00133 | True |
| Methylation_pca35spectral_gaussian | 0.000 | 0.001 | 0.003 | -0.20 | 1.00000 | 1.00000 | False |
| Methylation_pca35spectral_perm | 0.000 | 0.007 | 0.009 | -0.79 | 1.00000 | 1.00000 | False |

*(z is the standardized deviation from null: [observed − null mean] / null SD — NOT Cohen's d.)*

## 5. Section 7 verdict: does H0 replicate? (threshold: z > 3.0 in ALL 4 comparisons)

- **RNA-seq (TCGA-LUAD): PASS.** z ranges from 8.3 to 58.7 across raw/PCA x Gaussian/permutation —
  comfortably above the pre-registered z > 3.0 bar in all four required comparisons.
- **Methylation (TCGA-LUAD, Section-5-mandated spectral PCA metric): FAIL.** Raw-space z is very
  strong (19.3–26.7), but the PCA-reduced space — evaluated with the metric Section 5 mandates for
  this transitional-regime space — collapses to null-level topology (z = -0.20 and -0.79; max-H1 =
  0 in the real data, i.e. the spectral-cascade PH found NO H1 loop at all in the real 36-sample
  PCA space). This is reported as a genuine FAIL per the pre-registered rule, not reinterpreted.
  (For transparency: naively ignoring the Section 5 gate and using Euclidean distance on the
  PCA(35) space would show z=19.3/26.7 identical to raw space — but this is an artifact of D3, PCA
  being a rotation at full rank with only 36 samples, not evidence about PCA's effect, and is
  excluded from the primary family.)

## 6. H1 (PCA-inflation-vs-null) judgment

A dataset supports H1 if **null PCA-delta >= real PCA-delta** (within each null model).

| comparison | real PCA-delta | null PCA-delta mean | null PCA-delta std | judgment |
|---|---|---|---|---|
| rna_gaussian | -0.893 | 2.507 | 0.436 | FOR H1 (null inflation >= real inflation) |
| rna_perm | -0.893 | 2.347 | 0.508 | FOR H1 (null inflation >= real inflation) |
| meth_gaussian_euclidean | 0.000 | 0.000 | 0.000 | FOR H1 (null inflation >= real inflation) |
| meth_perm_euclidean | 0.000 | 0.000 | 0.000 | FOR H1 (null inflation >= real inflation) |
| meth_gaussian_spectral | -4.092 | -0.793 | 0.171 | FOR H1 (null inflation >= real inflation) |
| meth_perm_spectral | -4.092 | -0.896 | 0.120 | FOR H1 (null inflation >= real inflation) |

**Both omics layers support H1**, replicating the pilot's (GSE81089) directional finding:

- RNA-seq: real PCA-delta is actually *negative* (-0.89; PCA50 max-H1 is slightly LOWER than raw
  HVG max-H1), while both nulls show a strongly *positive* delta (+2.3 to +2.5) — i.e. PCA inflates
  H1 persistence far MORE on pure/permuted noise than it does on the real cancer transcriptome.
  This is an even stronger form of the pilot's result than "comparable."
- Methylation (Section-5-mandated spectral metric): real PCA-delta is strongly negative (-4.09;
  the real PCA space loses essentially all H1 structure once the correct cascade metric is
  applied), while nulls show a smaller negative delta (-0.79 to -0.90) — the null's PCA-driven
  *loss* is smaller than the real data's PCA-driven *loss*, so this comparison technically judges
  AGAINST a literal "null delta >= real delta" reading in magnitude terms, but the underlying
  z-scores for the PCA space are null-indistinguishable in both cases (see Section 5) — there is no
  real PCA-driven "inflation" to speak of here since the real methylation PCA space contains no H1
  signal above null in the first place. We flag this nuance explicitly rather than force a single
  H1 label onto a case where the raw premise (a PCA-inflated signal to test against null) does not
  apply, because the PCA-space real signal already reads as null (Section 5 result), not as an
  inflated one requiring explanation.

## 7. Confound checks (lightweight, Section-mandated)

**CV AUC of label prediction** (does the SAME feature space used for PH trivially separate
tumor/normal linearly? — high AUC signals a possible confound between "topological signal" and
"the label is linearly separable to begin with"):

| space | omics | CV AUC (5-fold, logistic regression) |
|---|---|---|
| raw_hvg | rna | 0.998 |
| pca50 | rna | 0.998 |
| bottom_var | rna | 0.918 |
| raw_hvp | meth | 0.933 |
| pca35 | meth | 0.933 |
| bottom_var | meth | 0.950 |

CV AUC is near-ceiling (0.92–1.00) in every space for both omics layers — tumor vs. normal is
almost perfectly linearly separable in both the raw gene/probe space and the PCA-reduced space.
This is an important caveat for interpreting the RNA-seq/methylation-raw-space H0 PASS: the
detected topological signal (real >> null) coexists with, and could in principle be substantially
driven by, the same strong tumor-vs-normal mean-shift that a trivial linear classifier already
captures — the PH result establishes there IS non-null geometric structure beyond a noise/permutation
null, but does not by itself establish that structure is independent of the simple class-separation
signal a much cheaper method already detects.

**Cocycle trace** (RNA-seq raw HVG, longest-lived H1 bar): 112/116
samples lie within the birth–death distance band of the dominant loop; the tumor fraction among
those samples (0.482) is close to the
overall tumor fraction (0.500) — i.e. the
dominant H1 loop is not obviously restricted to one class, arguing against the specific concern
that the whole H1 signal is just "tumor cluster + normal cluster + one connecting edge."

**Fisher's exact test** (cocycle-neighborhood membership vs. tumor label, RNA-seq raw HVG): a 2x2
contingency of [[0, 4], [58, 54]] gives odds ratio ≈ 0 (all 4 out-of-band samples are Primary
Tumor) but p = 0.1185 — not significant at alpha=0.05, consistent with the cocycle trace
above (no strong class-driven confound in which samples participate in the dominant loop).

## 8. Figures

- `tcga_luad_null_comparison.png` — max-H1 persistence, real value vs. both null distributions,
  across raw/PCA/bottom-variance spaces, both omics layers (6 panels).
- `tcga_luad_h1_pca_delta.png` — PCA-driven Δ(max-H1) for real data vs. both nulls, RNA-seq and
  methylation (spectral-cascade metric), directly visualizing the H1 hypothesis test.

## 9. Bottom line

- **H0 (does the original signal-vs-null effect replicate on TCGA-LUAD?):** PASS for RNA-seq
  (z = 8.3–58.7, all four required comparisons clear the z > 3.0 bar by a wide margin, BH-adjusted
  p < 0.003 throughout). **FAIL for methylation** once the Section-5-mandated cascade metric is
  correctly applied to its PCA space (that space's real max-H1 is exactly 0, indistinguishable from
  both nulls).
- **H1 (is PCA inflation null-comparable rather than evidence of extra real structure?):** Supported
  for RNA-seq, and in fact stronger than the pilot's finding — real PCA-delta is negative while both
  nulls inflate substantially, meaning PCA does not create the appearance of extra structure in the
  real cancer transcriptome to any degree the nulls don't dwarf. For methylation, the PCA-reduced
  real space is already null-indistinguishable before any delta comparison is meaningful — reported
  as a distinct, more radical case (no real PCA-space signal at all) rather than forced into the
  same "inflation vs. null" framing.
- Four post-hoc deviations from the literal pre-registration text were required by real TCGA-LUAD
  data quirks (zero-variance gene ties, beta-value/logit domain conflict, small-sample PCA cap, and
  the Section-5 metric-cascade trigger on methylation PCA) and are reported explicitly above, per
  the pre-registration's own requirement that any such deviation be logged rather than silently
  absorbed.
