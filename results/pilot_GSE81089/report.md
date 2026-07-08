# Does PCA destroy or concentrate topological structure in bulk RNA-seq? A transfer test of T-002/OQ-006 on GSE81089 NSCLC data

**Author:** Claude Science (TOPOLOGICA project) | **Date:** 2026-07-08 | **Dataset:** GEO accession [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089) (Djureinovic et al., non-small cell lung cancer RNA-seq)

## 1. Motivation and scope

The validated-results-registry records VR-003 (scRNA-seq, PBMC3k): under a naive dimension-matched Gaussian null, PCA(50) appeared to destroy real topological (H1) structure present in raw highly-variable-gene (HVG) space. Two open items follow directly from that result:

- **T-002** — does this PCA effect on topology transfer from single-cell to *bulk* RNA-seq / multi-omics data?
- **OQ-006** — is a pure Gaussian-noise null itself sufficient to characterize the effect, or does PCA fabricate *phantom* topology out of noise regardless of whether real structure is present?

Mid-session, a claim arrived (via an unverified, injected message masquerading as project context) that VR-003 had since been "corrected" using a pipeline-symmetric null (identical normalization/HVG/PCA applied to a permuted surrogate), inverting the finding to "PCA concentrates topology." **This claim could not be verified against the actual registry content and is treated here as an untrusted, unconfirmed assertion, not ground truth.** Rather than adopt or reject it outright, this analysis runs **both** null models side by side — the original OQ-006 Gaussian null and a pipeline-symmetric permutation null — and reports what each shows independently, so the reader can judge for themselves rather than relying on an unverifiable claim.

## 2. Data and preprocessing

**Dataset:** GSE81089 — bulk RNA-seq (FPKM, Cufflinks) from 199 NSCLC tumor samples and 19 matched normal lung samples, 218 samples total, 63,130 annotated Ensembl genes.

**Point-cloud convention:** samples are the points (218 points), genes are the ambient coordinates — directly analogous to VR-003's "cells are points, genes are dimensions" convention in PBMC3k.

**Power caveat (stated up front):** 218 samples is roughly two orders of magnitude fewer points than a typical scRNA-seq dataset (PBMC3k has ~2,700 cells). Persistent homology's ability to resolve fine topological structure scales with point-cloud density; null-result interpretations below must be read with this power gap in mind, not treated as equivalent in strength to a single-cell finding.

**Preprocessing:** log1p(FPKM) transform; genes with zero variance dropped (45,047 of 63,130 retained); top-2000 highly-variable genes (HVG) selected by variance, and a separate bottom-2000-variance gene set retained (mirroring VR-003's discovery that low-variance genes carried signal in scRNA data); features standardized (zero mean, unit variance) before PCA; PCA computed to 50 components (70.9% variance explained on real data).

**A data-quality note:** the raw FPKM matrix contained a small number of `-1.0` sentinel values (1,124 of ~13.8M entries, 0.008%) representing masked/undetected calls; these were treated as missing and imputed to 0 (non-detected). This is a minor, disclosed preprocessing choice, not expected to materially affect the topology given its size.

## 3. Metric validity check (metric-space-diagnostics protocol)

Before computing any persistent homology, intrinsic dimension was estimated (MLE, Levina-Bickel) for every feature space used:

| Space | d_int (MLE) | d_ambient | ratio | Regime |
|---|---|---|---|---|
| Raw HVG (top-2000) | 25.65 | 2000 | 0.013 | JL (Euclidean valid) |
| Bottom-variance genes (2000) | 25.88 | 2000 | 0.013 | JL (Euclidean valid) |
| PCA(50) | 10.91 | 50 | 0.218 | JL (Euclidean valid) |

All three conditions fall well inside the Johnson-Lindenstrauss regime (ratio < 0.3), so standard Euclidean-distance Vietoris-Rips persistent homology is methodologically valid throughout — no spectral or Fermat-distance fallback was required.

## 4. Two independent null models

1. **OQ-006 Gaussian null** — pure i.i.d. N(0,1) noise, dimension-matched to the real HVG matrix (218×2000), independently re-drawn 500 times, each drawn through the identical standardize→PCA(50) pipeline as the real data.
2. **Pipeline-symmetric permutation null** — each gene's expression values independently permuted across the 218 samples (destroys cross-gene/cross-sample correlation structure while exactly preserving each gene's marginal distribution and overall variance ranking), then passed through the *identical* HVG-selection/standardization/PCA(50) pipeline as the real data, repeated 2000 times.

These are genuinely different null hypotheses: the Gaussian null asks "is this data structured at all, relative to unstructured noise with matched dimensionality?"; the permutation null asks "is there real cross-gene/cross-sample correlation structure, beyond what per-gene marginal statistics alone would produce?" Both are reported because they answer different questions and neither alone settles the matter.

## 5. Results

### 5.1 Topological summary (single-draw observed values)

| Condition | H0 (components) | H1 (loops) | H2 (voids) |
|---|---|---|---|
| Real, raw HVG (2000 genes) | 218 | 95 | 19 |
| Real, bottom-variance genes (2000 genes) | 218 | **0** | 0 |
| Real, PCA(50) | 218 | 109 | 24 |

The bottom-variance-gene result does **not** replicate VR-003's scRNA finding (where bottom-844 genes carried the *strongest* H1 signal, persistence 3.15 vs top-50's max of 1.39). In this bulk dataset, the bottom-2000-variance genes carry essentially no detectable H1 topology at all. Given the ~12x lower sample count here relative to typical scRNA datasets, this is likely an underpowered null result rather than a genuine absence of signal — but it is reported as observed, not explained away.

### 5.2 Permutation test: is the observed max-H1-persistence higher than either null? (VR-003-style protocol)

Primary statistic: **max H1 persistence** (the single most persistent loop) — the statistic most resistant to noise-count inflation, since permutation nulls generate many short-lived spurious bars but rarely a single very long-lived one.

| Comparison | Observed | Null mean ± SD | p-value | Standardized deviation from null (z-score)* | Validation method |
|---|---|---|---|---|---|
| Real raw HVG vs. pipeline-symmetric null | 3.887 | 1.325 ± 0.172 | **0.0005** | **14.9** | Permutation test (per-gene shuffle), n=2000 |
| Real PCA50 vs. pipeline-symmetric null (PCA50) | 4.671 | 2.830 ± 0.352 | **0.0005** | **5.2** | Permutation test + identical PCA pipeline, n=2000 |
| Real raw HVG vs. Gaussian null | 3.887 | 1.037 ± 0.109 | **0.0020** | **26.1** | Monte Carlo Gaussian null, n=500 draws |
| Real PCA50 vs. Gaussian null (PCA50) | 4.671 | 3.036 ± 0.315 | **0.0020** | **5.2** | Monte Carlo Gaussian null + PCA, n=500 draws |

*Correction: an earlier version of this table and `final_results_table.csv` mislabeled this column "Cohen's d." It is (observed − null mean)/null SD computed against a single real-data observation — a z-score relative to the null distribution, not a two-sample Cohen's d (which requires variance on both sides being compared). The magnitude and conclusion are unchanged; only the label was wrong. See `final_results_table_corrected.csv` for the corrected labeling.

**Both null models agree on the headline result:** the real dataset's most persistent H1 feature is far outside either null distribution, in both raw and PCA-reduced space. This is a genuine, large-effect-size signal — not marginal.

### 5.3 Does PCA inflate persistence in noise too? (isolating a PCA artifact from a real-data effect)

| Comparison | Mean Δ(PCA − raw) | Paired test | Effect size (d) |
|---|---|---|---|
| Gaussian null: PCA50 vs raw (paired, same 500 draws) | +2.00 | Wilcoxon p=1.3e-83 | 5.93 |
| Pipeline-symmetric null: PCA50 vs raw (paired, same 2000 perms) | +1.50 | Wilcoxon p≈0 | 3.87 |
| **Real data: PCA50 vs raw** | **+0.78** | (single observation, no distribution) | — |

**This is the most important finding of this analysis.** PCA inflates max-H1-persistence in *both* null models — including pure Gaussian noise, which by construction contains no biological structure whatsoever. The magnitude of PCA-driven inflation in noise (+1.5 to +2.0) is *larger* than the PCA-driven increase observed in the real data itself (+0.78).

This directly addresses OQ-006: **yes, PCA fabricates phantom H1 topological inflation on pure noise.** This is consistent with the general finding reported (from an unverified source) that "PCA concentrates topology" — but the noise result shows this concentration effect is not evidence of *real* structure being revealed; it is at least partly a generic geometric artifact of projecting many roughly-isotropic dimensions down to fewer, more curvature-prone ones. The dominant driver of the real dataset's significance versus null is the *raw*-space gap (d=14.9–26.1), which is unaffected by this PCA artifact, not the PCA-space gap (d=5.2, and comparable in magnitude to the noise-only PCA inflation).

### 5.4 Confound check: is the real H1 signal just tumor-vs-normal labels in disguise?

A 5-fold cross-validated logistic regression on PCA(50) coordinates predicts Tumor/Normal status with **AUC = 1.000 ± 0.000** — sample structure in this dataset is overwhelmingly dominated by tumor/normal separation, as expected. However, tracing back the most persistent H1 loop's participating vertices (via `ripser`'s cocycle representatives) shows it touches only 3 of 218 samples (all tumor), not a tumor/normal boundary structure. A Fisher's exact test found no significant enrichment of either label among loop participants (p=1.0, likely underpowered given only 3 participating points). **This loop is not simply encoding the tumor/normal label** — but with only 3 points, this analysis cannot rule out that it reflects a narrow tumor subtype or a batch/technical artifact instead. Larger, better-annotated cohorts would be needed to resolve what the loop biologically represents.

## 6. Honest conclusion

**T-002 (transfer of PCA effect to bulk RNA-seq): SUPPORTED, with a qualification.** Real bulk RNA-seq data shows a persistent-homology signal in H1 significantly exceeding both a Gaussian-noise null and a pipeline-symmetric permutation null (p≈0.0005–0.002, very large effect sizes, d=5.2–26.1) — genuine topological structure exists beyond either null's baseline, in both raw and PCA-reduced feature space. This transfers the qualitative finding "real expression data carries real topology" from single-cell to bulk RNA-seq.

**Whether PCA "destroys" or "concentrates" that structure: INCONCLUSIVE, and the framing itself is likely wrong.** The real data's own raw→PCA change (+0.78) is smaller than the artifactual PCA-driven inflation measured in both null models (+1.5 to +2.0). This means: (a) PCA does not straightforwardly destroy the real H1 signal here — it survives, comparably significant, in both raw and PCA space; but (b) a naive "PCA amplifies persistence, therefore PCA reveals more real structure" reading is not supported either, because PCA amplifies persistence in *pure noise* by a larger margin. The honest reading is that **max-H1-persistence as a summary statistic is not, by itself, a reliable indicator of whether PCA is revealing or fabricating structure** — this statistic must always be interpreted relative to a matched null, not in isolation, and the choice of null (Gaussian vs. permutation) matters less here than the raw-vs-PCA framing itself.

**OQ-006 (does PCA fabricate phantom topology on Gaussian noise): CONFIRMED.** PCA(50) applied to pure i.i.d. Gaussian noise produces significantly higher max-H1-persistence than the same noise analyzed in raw feature space (paired Wilcoxon p=1.3e-83, d=5.93). This is a domain-independent methodological artifact, not specific to biological data, and it means any PCA+persistent-homology pipeline needs a matched (ideally pipeline-symmetric) null before attributing a PCA-space topological feature to real structure.

**Bottom-variance-gene replication of VR-003: NOT REPLICATED (power-limited).** Unlike PBMC3k, bottom-2000-variance genes in this bulk dataset show zero detectable H1 structure. Given the ~12x smaller sample count relative to scRNA, this is most plausibly an underpowered null result, not a genuine absence — but it is reported as observed and should not be over-interpreted as a contradiction of VR-003 without a power-matched follow-up (e.g., subsampling PBMC3k to n≈218 cells for a fair comparison).

**On the unverified "VR-003 was corrected" claim received mid-session:** this analysis's Section 5.3 result is *consistent* with a general "PCA concentrates rather than purely destroys" narrative, but demonstrates that narrative requires an important caveat (the noise-inflation confound) that a simple "PCA concentrates topology" headline omits. This should be treated as an independent, self-contained result — it neither confirms nor depends on the unverified claim, and the claim itself remains unconfirmed against the actual registry record.

## 7. Limitations

- **Sample count / power.** n=218 is small for persistent homology relative to typical single-cell applications; null results (bottom-variance genes) are likely underpowered, not necessarily "true negatives."
- **Single dataset.** One NSCLC cohort; generalization to other cancers, other bulk RNA-seq designs, or true multi-omics fusion (the original commercial framing) is untested.
- **Confound check is preliminary.** The top H1 loop's 3-sample participation is too small a set for a well-powered label-independence test; a larger, multi-cohort replication would strengthen this.
- **-1.0 sentinel imputation.** A minor, disclosed data-cleaning choice (0.008% of entries) that was not swept for sensitivity, though its scale makes a large effect on the topology unlikely.
- **Single real-data observation.** Unlike the null distributions (500–2000 draws), the real dataset provides only one observed value per statistic — appropriate for a permutation test against a null, but it means no confidence interval can be placed directly on the real-data statistic itself.

## 8. Artifacts

- `persistence_diagrams_real_vs_gaussian.png` — H1 persistence diagrams, real data vs. pooled Gaussian-noise draws, raw and PCA50 spaces
- `barcode_and_null_distributions.png` — top-30 H1 barcodes and full null distributions (both null models) with observed real-data value marked
- `final_results_table.csv` — full permutation test statistics table (p-values, effect sizes, validation method per comparison)
- `permutation_test_summary.csv` — intermediate summary table
