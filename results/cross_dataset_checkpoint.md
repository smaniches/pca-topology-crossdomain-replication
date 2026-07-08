# Cross-Dataset Replication Checkpoint

**Status:** Post-Phase-2 (replication complete for 3 of 3 planned confirmatory datasets), pre-ablation.
**Pre-registration:** locked 2026-07-08T13:50:55Z, SHA-256 `5e539309747188a2e77aa36bf1f9aecd3b50b8493803efa756ee4715773f0517`
**Purpose of this document:** record, per dataset, fetch/analysis status, BH-adjusted p-values, z-scores,
and effect-size-threshold pass/fail, exactly as computed — before any ablation or confound-audit work
touches these results. This is the "checkpoint before Phase 2 touches new data" record the pre-registration
required; it is written after the 3 confirmatory fetches/analyses complete and before ablation begins.

## 1. Datasets analyzed

| # | Dataset | Modality | Cancer type / tissue | N samples | Fetch status | Analysis status |
|---|---|---|---|---|---|---|
| Prior (excluded from confirmatory family) | GSE81089 | Bulk RNA-seq | NSCLC (lung) | 218 (199T/19N) | done (pilot, pre-dates pre-registration) | done |
| 1 | GSE146889 | Bulk RNA-seq | Colorectal/endometrial/ovarian (MMR-deficiency) | 176 (91T/85N) | COMPLETE | COMPLETE |
| 2 | CPTAC-CCRCC (PDC000127) | Proteomics (protein abundance) | Renal | 194 (110T/84N) | COMPLETE, no substitution needed | COMPLETE |
| 3 | TCGA-LUAD (GDC) | RNA-seq + DNA methylation | Lung adenocarcinoma | 116 (RNA-seq); 36 (paired methylation subset) | COMPLETE via GDC API | COMPLETE, 4 declared deviations |

All three planned confirmatory datasets were fetched successfully; none required a substitution or were
blocked. No synthetic or simulated data was used at any point — a fetch or analysis failure would have
been reported as BLOCKED per the sub-agent task instructions, not backfilled with fabricated values.

## 2. Full confirmatory test family (BH-FDR correction, alpha=0.05)

Per pre-registration Section 6, the correction spans the primary raw-space/PCA-space x
Gaussian-null/permutation-null family across all confirmatory datasets and omics layers (16 tests total —
the pilot's own tests are excluded, consistent with its role as the prior that motivated pre-registration,
not a confirmatory test of itself).

| dataset | layer | condition | z | p (raw) | p (BH) | reject@0.05 | z>3.0 | H0 verdict |
|---|---|---|---|---|---|---|---|---|
| GSE146889 (2nd RNA-seq cohort) | RNA-seq | Real raw HVG | 44.87 | 0.00200 | 0.00246 | True | True | PASS |
| GSE146889 (2nd RNA-seq cohort) | RNA-seq | Real raw HVG | 39.59 | 0.00050 | 0.00114 | True | True | PASS |
| GSE146889 (2nd RNA-seq cohort) | RNA-seq | Real PCA50 | 12.02 | 0.00200 | 0.00246 | True | True | PASS |
| GSE146889 (2nd RNA-seq cohort) | RNA-seq | Real PCA50 | 11.37 | 0.00050 | 0.00114 | True | True | PASS |
| CPTAC-CCRCC (proteomics) | Protein abundance | raw HVG vs Gaussian null | 23.18 | 0.00200 | 0.00246 | True | True | PASS |
| CPTAC-CCRCC (proteomics) | Protein abundance | PCA50 vs Gaussian null | 4.40 | 0.00399 | 0.00456 | True | True | PASS |
| CPTAC-CCRCC (proteomics) | Protein abundance | raw HVG vs pipeline-symmetric null | 21.39 | 0.00050 | 0.00114 | True | True | PASS |
| CPTAC-CCRCC (proteomics) | Protein abundance | PCA50 vs pipeline-symmetric null | 4.64 | 0.00050 | 0.00114 | True | True | PASS |
| TCGA-LUAD | RNA-seq | RNAseq_raw_gaussian | 58.67 | 0.00200 | 0.00246 | True | True | PASS |
| TCGA-LUAD | RNA-seq | RNAseq_raw_perm | 47.53 | 0.00050 | 0.00114 | True | True | PASS |
| TCGA-LUAD | RNA-seq | RNAseq_pca50_gaussian | 9.19 | 0.00200 | 0.00246 | True | True | PASS |
| TCGA-LUAD | RNA-seq | RNAseq_pca50_perm | 8.34 | 0.00050 | 0.00114 | True | True | PASS |
| TCGA-LUAD | Methylation | Methylation_raw_gaussian | 19.34 | 0.00200 | 0.00246 | True | True | PASS |
| TCGA-LUAD | Methylation | Methylation_raw_perm | 26.69 | 0.00050 | 0.00114 | True | True | PASS |
| TCGA-LUAD | Methylation | Methylation_pca35spectral_gaussian | -0.20 | 1.00000 | 1.00000 | False | False | FAIL |
| TCGA-LUAD | Methylation | Methylation_pca35spectral_perm | -0.79 | 1.00000 | 1.00000 | False | False | FAIL |

14 of 16 tests PASS at BH-adjusted alpha=0.05 and clear the pre-registered z>3.0 effect-size bar. The 2
FAILs are both TCGA-LUAD methylation PCA-space comparisons (Section-5-mandated spectral metric), where the
real data's PCA-reduced methylation space contains literally zero H1 structure (max-H1 = 0), indistinguishable
from both nulls.

## 3. Per-dataset H0 verdict (does "real signal exceeds null" replicate?)

| Dataset | Layer | H0 verdict | Notes |
|---|---|---|---|
| GSE146889 | RNA-seq | **PASS** | z = 11.4–44.9 across all 4 required comparisons |
| CPTAC-CCRCC | Proteomics | **PASS** | z = 4.4–23.2 across all 4 required comparisons |
| TCGA-LUAD | RNA-seq | **PASS** | z = 8.3–58.7 across all 4 required comparisons |
| TCGA-LUAD | Methylation | **PARTIAL — raw-space PASS, PCA-space FAIL** | Raw-space z = 19.3–26.7 (strong); PCA-space (spectral-cascade metric, Section-5-mandated) collapses to null-indistinguishable (z = -0.20, -0.79), real max-H1 = 0 |

**3 of 4 layers show full H0 replication. The 4th (TCGA-LUAD methylation) shows partial replication: strong
signal in raw feature space, but no detectable H1 structure survives PCA-reduction once the correct
(spectral-cascade) distance metric is applied, as mandated by the intrinsic-dimension gate for that
36-sample, transitional-regime space.** This is reported as a genuine, non-overridden divergence from the
other three layers, not smoothed into a uniform "replicates everywhere" narrative.

## 4. Per-dataset H1 verdict (is PCA-inflation null-comparable, not evidence of extra real structure?)

| Dataset | Layer | H1 verdict | Real PCA-delta | Null PCA-delta range |
|---|---|---|---|---|
| Pilot (GSE81089) | RNA-seq | Supported | +0.78 | +1.50 to +2.00 |
| GSE146889 | RNA-seq | Supported | +1.36 | +2.08 to +2.18 |
| CPTAC-CCRCC | Proteomics | Supported | +0.85 | +2.01 to +2.10 |
| TCGA-LUAD | RNA-seq | Supported (stronger form) | **−0.89** (negative) | +2.35 to +2.51 |
| TCGA-LUAD | Methylation | Supported, but with a distinct interpretation (see below) | −4.09 | −0.90 to −0.79 |

**H1 replicates in every layer tested, across four independent real-data sources and two omics modalities
(transcript and protein abundance), plus a partial fifth (methylation) with a different mechanism.** In
every case, the null models' PCA-driven persistence change exceeds (or in TCGA-LUAD RNA-seq, dramatically
exceeds — real delta is negative) the real data's own PCA-driven change. TCGA-LUAD methylation is a
qualitatively different case: the real PCA-space signal is already null-indistinguishable before any
"inflation" framing applies (Section 5's own H0 result establishes there's no real PCA-space structure to
compare against a null in the first place), so its H1 entry is reported with that caveat rather than folded
into the same "inflation exceeds null" statistic as the other four.

## 5. Confound-check summary (lightweight, pre-registered)

| Dataset | CV AUC (label separability) | Dominant-loop label association |
|---|---|---|
| GSE81089 (pilot) | 1.000 | Not significant (Fisher p=1.0, underpowered, 3 participants) |
| GSE146889 | not reported here — see full report | Significant (Fisher p=2.9e-09, 45/176 samples) |
| CPTAC-CCRCC | not reported here — see full report | Not significant (Fisher p>0.13, underpowered, 3-5 participants) |
| TCGA-LUAD | 0.918–0.998 (near-ceiling in every space) | Not significant (Fisher p=0.12) |

Tumor/normal class separability is near-ceiling in every dataset tested (CV AUC 0.92–1.00), meaning the
detected topological signal coexists with a very strong, simple linear mean-shift between classes. The
PH-based tests establish there IS non-null geometric structure beyond a matched permutation/noise null, but
do NOT by themselves establish that this structure is independent of the same class-mean-shift signal a
much cheaper linear classifier already captures. This coexistence is flagged consistently across datasets
as a caveat for the confound-attribution-audit phase to address directly (LOCO, incremental decomposition,
residualization against the label).

## 6. Post-hoc deviations disclosed

Per the pre-registration's own requirement, any necessary deviation from the literal locked protocol is
logged here rather than silently absorbed:

- **TCGA-LUAD D1:** Bottom-2000-variance genes selected among genes with variance > 0 (5,556/60,660 genes
  have exactly zero variance in this cohort; the literal bottom-2000 would be a degenerate all-identical
  point cloud).
- **TCGA-LUAD D2:** Methylation missing-value handling deviated from the literal "missing → 0" rule
  (incompatible with logit(beta), which diverges at 0); values clipped to [1e-6, 1-1e-6] before the logit
  transform.
- **TCGA-LUAD D3:** Methylation PCA capped at 35 components (not the pre-registered 50), because only 36
  samples are available in the matched methylation subset (sklearn's full-SVD rank cap). This makes
  methylation's Euclidean-PCA space a rotation-only, full-rank isometry of the raw space — not an
  independent test — and is reported for transparency only, excluded from the primary family.
- **TCGA-LUAD D4:** Per Section 5's intrinsic-dimension gate, methylation PCA space (d_int/d_amb ratio =
  0.320, crossing the 0.3 transitional threshold) required the spectral-cascade metric instead of raw
  Euclidean distance; this is the metric used for the primary-family methylation-PCA comparisons.

No deviations were required for GSE146889 or CPTAC-CCRCC; both datasets were analyzed under the literal
pre-registered protocol without modification.

## 7. What this checkpoint does and does not establish

**Established with reasonable confidence:**
- Real, disparate omics data (bulk RNA-seq from 3 independent cohorts across 4 cancer types/tissue groups,
  plus one independent proteomics cohort) carries H1 persistent-homology structure that significantly
  exceeds both a Gaussian-noise null and a pipeline-symmetric permutation null, at a pre-registered
  z>3.0 threshold, in raw feature space essentially universally, and in PCA-reduced space in the large
  majority of layers tested.
- The pilot's central methodological caution — that PCA-driven increases in max-H1-persistence are, by
  themselves, an unreliable signal of "PCA revealing real structure" because null models show comparable or
  larger increases — replicates across every layer tested, including a stronger form (negative real delta)
  in TCGA-LUAD RNA-seq.

**Not yet established, deferred to later phases:**
- Whether the detected topological signal is independent of the tumor/normal class mean-shift (near-ceiling
  CV AUC in every dataset) — this requires the confound-attribution-audit's LOCO/incremental-decomposition/
  residualization protocol, not yet run.
- Whether the specific hyperparameters (2000 HVG, 50 PCs, imputation-to-zero) materially drive these
  results — this requires the ablation sweep, not yet run.
- Robustness of the TCGA-LUAD methylation-PCA divergence (n=36 is small; a larger paired-methylation cohort
  would strengthen or weaken this specific finding).

## 8. Next steps

1. Full ablation sweep (gene count, PC count, imputation rule, permutation count) per
   universal-ablation-engine, on at least the pilot + one confirmatory dataset.
2. Full confound-attribution-audit (LOCO, incremental decomposition, hard-negative/within-stratum/
   residualization controls, block-bootstrap CIs) to directly address the CV-AUC/topology coexistence
   flagged in Section 5 above.
3. Synthesis of all of the above into a final, honest verdict for the paper.
