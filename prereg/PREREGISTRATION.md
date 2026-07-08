# PRE-REGISTRATION DOCUMENT
## Cross-domain replication test: does PCA inflate persistent-homology signal beyond a matched null?

**Status:** LOCKED — no changes permitted after this document is timestamped and checksummed below.
**Scope:** This document governs all NEW dataset analyses in this program. The GSE81089 pilot
(already analyzed, reported in report.md) is NOT covered by this pre-registration — it was
exploratory/pilot work that motivated this locked protocol, and is treated as a prior, not as
one of the pre-registered confirmatory tests below.

---

## 1. Primary hypothesis (stated directionally, per the pilot's actual finding)

**H1 (primary, confirmatory):** In each new dataset, PCA(50) inflates max-H1-persistence
relative to BOTH null models (i.i.d. Gaussian noise, and pipeline-symmetric per-gene/per-feature
permutation) by an amount that is COMPARABLE TO OR LARGER THAN the PCA-driven increase observed
in the real data itself.

This is the pilot's actual finding (GSE81089: real PCA-delta = +0.78, Gaussian-null PCA-delta =
+2.00, pipeline-null PCA-delta = +1.50 — the null inflation exceeded the real inflation). The
pre-registered claim is therefore: **"PCA-driven persistence increase is not, by itself, evidence
of real structure being revealed — because a matched null shows the same or larger increase from
pure or permuted noise."** A dataset is judged FOR this hypothesis if null PCA-delta >= real
PCA-delta (within each null model separately); AGAINST it if real PCA-delta clearly exceeds both
null PCA-deltas.

**H0 (secondary, does the ORIGINAL signal-vs-null effect replicate):** Real data max-H1-persistence
significantly exceeds the null distribution (both Gaussian and pipeline-symmetric), in both raw
and PCA-reduced feature space. This is judged against the effect-size threshold in Section 5.

These are two distinct, separately-reported claims. A dataset can support H0 (real structure
exists beyond null) while also supporting H1 (PCA inflation is a null-comparable artifact, not
evidence of MORE real structure being revealed by PCA specifically).

## 2. Primary statistic

**Max H1 persistence**: the birth-death gap of the single longest-lived H1 (1-dimensional loop)
feature in the persistent homology barcode. This was chosen (in the pilot, retroactively — see
scope note above) because it is the statistic most resistant to bar-count inflation from short-
lived permutation-null noise.

## 3. Fixed hyperparameters (identical across all datasets unless Section 4 modality exception applies)

- HVG (highly-variable-gene) count: **2000**
- PCA components: **50**
- Permutation count, pipeline-symmetric null: **2000** iterations
- Monte Carlo draw count, Gaussian null: **500** draws
- Imputation rule: missing / sentinel values -> **0** (treated as non-detected)
- Random seed: **42** (all stochastic steps)
- Standardization: zero-mean, unit-variance scaling applied before PCA, fit on the condition being tested (no test-fold leakage across null draws)

## 4. Per-modality normalization rules (locked now, before any modality-specific data is fetched)

- **Bulk RNA-seq** (2nd GEO cohort, TCGA RNA-seq): log1p(FPKM or TPM), matching the pilot.
- **TCGA second omics layer** (methylation or CNV, whichever is more readily accessible via GDC):
  methylation -> logit(beta-value); CNV -> log2(copy-ratio). Then identical HVG-selection /
  standardization / PCA / PH steps as RNA-seq.
- **CPTAC proteomics** (intensity data): median-normalize across samples, then log2 transform.
  Then identical HVG-selection / standardization / PCA / PH steps.
- **No other pipeline step may vary by dataset or modality.** Any additional deviation required
  by a real data quirk (e.g. a new missingness pattern) must be reported explicitly as a
  post-hoc deviation, not silently absorbed into "the same pipeline."

## 5. Intrinsic-dimension gate (carried into every new dataset)

Before computing persistent homology on ANY new dataset/condition, estimate intrinsic dimension
via MLE (Levina-Bickel) and TwoNN (Facco et al.), per the metric-space-diagnostics protocol.
If d_int/d_amb >= 0.3 (outside the JL regime) for any condition, the Euclidean-Spectral-Fermat
metric cascade must be applied instead of plain Euclidean VR persistent homology, and this
deviation must be reported explicitly (which metric was used, why) rather than silently computed
with Euclidean distance regardless.

## 6. Multiple-testing correction

Benjamini-Hochberg (BH) FDR correction at alpha=0.05, applied across the full family of primary-
statistic tests: up to 5 datasets (pilot + up to 4 new) x 2 null models (Gaussian, pipeline-
symmetric) x 2 feature spaces (raw, PCA50) = up to 20 tests. Both raw and BH-adjusted p-values
are reported for every comparison. The pilot's own tests are included in this family for full
transparency, even though the pilot itself was exploratory (see scope note).

## 7. Effect-size threshold (practical significance, not just p < alpha)

A dataset is judged to REPLICATE H0 (real structure exceeds null) only if the standardized
deviation from null (z-score: [observed − null mean] / null SD) exceeds **3.0** in BOTH raw and
PCA-reduced space, against BOTH null models. This threshold was chosen to clear the pilot's
weakest single margin (PCA vs. pipeline-null z = 5.2) with room to spare, while remaining well
above conventional large-effect thresholds (Cohen's benchmark of d=0.8 for "large"). A dataset
not meeting this bar in all four required comparisons is reported as NOT replicating H0, without
reinterpretation after the fact.

## 8. Planned datasets

1. GSE81089 (pilot, already analyzed — reference prior, not a confirmatory test)
2. One additional, independent bulk RNA-seq cohort from GEO (distinct cancer/tissue type from GSE81089)
3. TCGA RNA-seq expression (via GDC API), plus a second omics layer (methylation or CNV) if accessible
4. CPTAC proteomics (via its public data portal), or a substitute real proteomics/multi-omics
   dataset (e.g. PRIDE, GEO proteomics series) if CPTAC access is blocked — substitution must be
   reported explicitly as a deviation from the original plan, never silently swapped in.

## 9. Decision rule for "does this replicate across datasets"

The overall program verdict is a simple count: report, for each dataset, PASS/FAIL against H0
(Section 7 threshold) and the H1 judgment (Section 1), with BH-adjusted p-values (Section 6).
No aggregate "replicated in X of Y" claim is treated as more than a descriptive summary — each
dataset's result stands on its own with its own reported statistics. Any dataset that cannot be
fetched or analyzed (e.g. blocked network access) is reported as BLOCKED, not silently dropped
from the denominator.

---

## Lock

This document is considered LOCKED at the timestamp and checksum below. No hyperparameter,
threshold, or hypothesis wording in Sections 1–9 may be altered after this point without an
explicit, separately-logged deviation note (which itself would compromise the confirmatory
status of any results obtained after the deviation).

**LOCKED AT (UTC):** 2026-07-08T13:50:55.244926Z
**SHA-256 CHECKSUM OF ABOVE CONTENT:** 5e539309747188a2e77aa36bf1f9aecd3b50b8493803efa756ee4715773f0517
