# Confound-Attribution Audit: Is the GSE146889 H1 Signal a Tumor/Normal Class-Mean-Shift Artifact?

**Project:** TOPOLOGICA PCA-topology-cross-domain replication | **Dataset:** GSE146889 (colorectal/endometrial/ovarian mismatch-repair-deficiency bulk RNA-seq, GEO)
**Protocol:** `confound-attribution-audit` skill, applied identically to the completed GSE81089 pilot audit
**Data refetched fresh** from `ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/suppl/GSE146889_GeneCount.tsv.gz` and reprocessed per the locked pre-registration's fixed hyperparameters (log1p(RPKM), drop-footer-row, drop zero-variance genes, HVG=2000, PCA=50, seed=42, standardize-then-PCA). Reproduced sample counts (91 tumor / 85 normal / 176 total), gene counts (64,253 after dropping the footer row; 51,634 after dropping zero-variance genes), PCA(50) variance explained (88.7%), real max-H1 persistence on PCA(50) (7.564), and 5-fold CV tumor/normal AUC (0.925 ± 0.019) all match `report_GSE146889.md` exactly, confirming the refetch/reprocess is faithful to the locked pipeline.

## 1. Named confound

**Tumor/normal class-mean-shift magnitude** — identical operationalization to the GSE81089 audit: the unit vector `u = (mean_tumor − mean_normal) / ‖mean_tumor − mean_normal‖` in standardized HVG(2000) space, and each sample's scalar projection onto `u`.

This dataset's own replication report already flagged a difference from the pilot before this audit began: the top H1 loop touches 45/176 samples and is significantly enriched for tumor status (Fisher exact p=2.87e-09) — unlike the pilot, where the loop touched only 3 samples with no label enrichment. This audit was run expecting the confound-independence result might not replicate, and reports the evidence as found.

## 2. Confound spectrum

| Quantity | Correlation with confound direction | Tier |
|---|---|---|
| PC1 | r = −0.889, p = 4.6e-61 | correlated |
| PC2 | r = 0.450, p = 3.7e-10 | correlated |
| PC3–PC50 | \|r\| < 0.06, p > 0.49 (all) | free |
| 1-D confound-direction alone | — | AUC = 0.973 using only this single projection |
| Variance explained by the confound direction | 22.3% of total HVG variance (446.6 / 2000.0) | — |

**Reading:** as in the pilot, the confound concentrates in PC1–PC2; PC3 onward are free of it. Here the confound explains a larger share of total variance (22.3% vs. the pilot's 11.2%) and PC1 alone correlates more strongly (r=−0.89 vs. r=0.87) — consistent with this dataset's somewhat stronger classifier separability signal relative to random noise, though its absolute CV AUC (0.925) is lower than the pilot's (1.000).

## 3. Frozen-baseline / group decomposition: within-class subsets vs. mixed set

| Condition | n | Observed max-H1 | z (Gaussian null) | p (Gaussian) | z (permutation null) | p (permutation) |
|---|---|---|---|---|---|---|
| Mixed (tumor+normal) | 176 | 7.564 | 12.02 | 0.0020 | 11.37 | 0.0005 |
| **Tumor-only** | 91 | **5.880** | **4.60** | 0.0033 | **4.24** | 0.0033 |
| **Normal-only** | 85 | **6.402** | **5.35** | 0.0033 | **5.09** | 0.0033 |

![Within-class null comparison]({artifact:within_class_null_comparison_GSE146889})

**This does not replicate the pilot's decisive first result.** In GSE81089, the tumor-only subset reproduced the mixed-set observed statistic *exactly* (3.887 = 3.887) with only a mild reduction in z (23.0 vs 26.1). Here, **neither single-class subset reproduces the mixed-set value**: tumor-only (5.880) is 78% of the mixed value (7.564), and normal-only (6.402) is 85% of it. Both single-class subsets still individually clear the pre-registered z>3.0 effect-size bar against both null models (z=4.2–5.3), so each class *does* carry detectable topological structure beyond a matched null — but the effect is markedly weaker than in the mixed set (z=4.6–5.3 vs. z=11.4–12.0), and weaker in relative terms than the pilot's near-parity between tumor-only and mixed (z=23.0 vs 26.1, a 12% reduction, versus a 62% reduction here for tumor-only). The mixed-set signal is not simply "the tumor-only signal, unchanged" — some of it depends on the two classes coexisting in the same point cloud.

## 4. Explicit control 1 — within-stratum effect (confound-quartile binning)

All 176 samples (tumor + normal) were binned into quartiles of the confound projection (distance along the tumor/normal mean-shift direction), and max-H1 persistence recomputed within each quartile against its own size-matched (n=44) Gaussian and pipeline-symmetric permutation nulls.

| Quartile | n | n_tumor | n_normal | Observed max-H1 | z (Gaussian) | z (permutation) |
|---|---|---|---|---|---|---|
| Q1 (least tumor-like) | 44 | 1 | 43 | 1.985 | 7.41 | 6.07 |
| Q2 | 44 | 9 | 35 | 2.722 | 12.16 | 9.80 |
| Q3 | 44 | 37 | 7 | 2.378 | 9.94 | 8.75 |
| Q4 (most tumor-like) | 44 | 44 | 0 | 4.491 | 23.59 | 19.55 |

![Within-stratum and residualization]({artifact:within_stratum_and_residualization_GSE146889})

**Every quartile clears the pre-registration's z>3.0 bar**, by a wide margin (min z=6.07, all ≥6.07). This part of the pilot's pattern **does replicate**: fixing the confound at a narrow band does not erase the topological excess over null. Note a structural difference from the pilot's stratification, however — because tumor/normal is highly (though not perfectly) separable along this exact projection, the quartiles are not label-balanced (Q1 is 98% normal, Q4 is 100% tumor); this control here is closer to "does the signal persist across the tumor-to-normal spectrum" than a clean within-class stratification, and Q4's strong z (23.6) is not fully separable from a pure tumor-class effect given Q4 contains zero normal samples.

## 5. Explicit control 2 — residualization on the confound direction

The linear class-mean-shift direction was regressed out exactly (each sample's class-specific mean shifted to the pooled grand mean, per-gene, in standardized HVG space), then the identical PCA(50) pipeline re-fit on the residualized matrix.

| Space | Observed max-H1 | z (Gaussian) | z (permutation) | 5-fold CV AUC (tumor/normal) | Top-loop composition | Fisher p (tumor enrichment) |
|---|---|---|---|---|---|---|
| Full HVG→PCA50 (class-mean intact) | 7.564 | 12.02 | 11.37 | 0.925 ± 0.019 | 45 samples (40 tumor / 5 normal) | 2.87e-09 |
| **Class-mean-residualized HVG→PCA50** | **8.415** | **14.35** | **13.52** | **0.206 ± 0.080** | 22 samples (10 tumor / 12 normal) | 0.650 |

![Bootstrap CI and AUC check]({artifact:bootstrap_ci_and_auc_check_GSE146889})

**Verification the residualization worked as intended:** class means are matched to numerical precision (max abs residual mean difference = 8.1e-16). 5-fold CV AUC collapses from 0.925 to 0.206 (a permutation test against a label-shuffled null, n=200, gives p=0.005 — the collapse itself is non-trivial, echoing the pilot's own residual-AUC finding, rather than evidence the residualization failed).

**The aggregate max-H1 signal survives residualization and even strengthens** (z=14.35 post- vs 12.02 pre-residualization, both against the identical Gaussian null) — this part of the pilot's pattern replicates: the topological excess over null is not annihilated by removing the exact linear quantity a classifier uses to separate the classes.

**However, the specific dominant loop's identity changes materially.** Pre-residualization, the top H1 loop touches 45 samples, 89% tumor, significantly enriched for tumor status (Fisher p=2.87e-09) — this is the same finding already flagged in `report_GSE146889.md`. Post-residualization, the new top loop touches a different, smaller set of 22 samples at a near-baseline 45%/55% tumor/normal split (Fisher p=0.650, no enrichment). This is a genuine departure from the pilot, where the single dominant loop's identity and composition were not examined pre/post-residualization because the pilot's loop showed no enrichment to begin with. Here, the *specific feature an analyst would point to as "the loop"* changes identity when the confound is removed — the aggregate statistic (max persistence, a single scalar) survives, but the geometric object generating it does not.

## 6. Block-bootstrap CI

Bootstrap resampled each condition's own sample set (mixed n=176, tumor-only n=91, normal-only n=85) with replacement, 2000 resamples, recomputing observed max-H1 persistence each draw and comparing against the fixed matched Gaussian null from Section 3/4 (held fixed across bootstrap draws, per the same efficiency approximation used in the pilot).

| Subset | n | Δ (point estimate) | 95% CI | Excludes 0? |
|---|---|---|---|---|
| Mixed (tumor+normal) | 176 | 4.380 | [2.547, 9.175] | **Yes** |
| **Tumor-only** | 91 | 2.294 | **[−0.695, 3.975]** | **No** |
| **Normal-only** | 85 | 2.781 | **[−1.645, 4.123]** | **No** |

**This is the most decisive divergence from the pilot.** In GSE81089, the tumor-only block-bootstrap CI excluded zero comfortably ([1.37, 3.77]). Here, **neither single-class subset's bootstrap CI excludes zero** — both the tumor-only and normal-only signal-beyond-null estimates are statistically indistinguishable from zero once the resampling variability of the observed statistic itself is accounted for, even though the same subsets individually cleared a p<0.0033 threshold in the fixed-null permutation test (Section 3). This is not a contradiction: the permutation test in Section 3 asks "is this one observed value unusual relative to a null distribution," which it is; the bootstrap in this section asks "if I resampled this same population of 91 (or 85) samples again, would I reliably see an excess over null," and the answer here is no — the point estimate is fragile to which particular samples are up- or down-weighted in the resample. The single-class signal in this dataset is real in the single observed draw but not robust in the stability sense that mattered most for the pilot's genuine-signal verdict.

## 7. Verdict

Applying the skill's reading rubric to each control separately, since this dataset does not give a uniform answer across all four:

**PARTIALLY GENUINE, WEAKER AND MORE CONFOUND-ENTANGLED THAN THE PILOT — this is a materially different, more qualified result than GSE81089's clean confound-independence finding, and should be reported as such rather than smoothed over.**

- **What replicates:** the within-stratum control (Section 4, all 4 quartiles z≥6.1, clearing the pre-registered z>3.0 bar) and the aggregate residualization control (Section 5, max-H1 signal *increases* in z after linear confound removal, from 12.0 to 14.3, while CV AUC collapses to nonspecific levels) both replicate the pilot's pattern. Taken alone, these two controls would support a "genuine, confound-independent" reading.
- **What does NOT replicate:** the within-class decomposition (Section 3) and the block-bootstrap CI (Section 6) — the two controls that gave the pilot's most decisive, least assumption-laden evidence — both fail to replicate here. Neither the tumor-only nor the normal-only subset reproduces the mixed-set statistic (unlike the pilot's exact match), and neither subset's bootstrap CI on the signal-beyond-null excludes zero (unlike the pilot's clean exclusion). The signal in each single-class subset is present in the one observed draw (Section 3's permutation test) but not stable under resampling (Section 6's bootstrap).
- **What is new and additional evidence of confound entanglement:** the dominant H1 loop in the intact (non-residualized) space is significantly enriched for tumor status (Fisher p=2.87e-09, already flagged in the underlying replication report) — the opposite of the pilot's finding of no enrichment. After residualization this specific enrichment disappears, but the loop's *identity* also changes — the aggregate scalar statistic surviving residualization does not mean the same geometric feature survives; a different loop, with different sample composition, takes over as the most persistent one.
- **Scope qualifier:** this dataset's H1 signal is best described as **partially attributable to the tumor/normal axis at the level of the single most-persistent feature, while the aggregate excess-over-null persists after linear confound removal.** It is not correct to say this dataset "confirms" the pilot's confound-independence finding, nor is it correct to say the signal is "purely" the confound (the residualized aggregate signal, and the within-stratum result, both argue against a pure-confound explanation). The honest statement is: **confound-independence is dataset-dependent** — strong and multiply-confirmed in GSE81089, present but substantially weaker and only partially confirmed in GSE146889, with the two controls that most directly test single-class robustness (within-class exact match, bootstrap CI) failing to replicate.

## 8. Limitations

- **Mixed tissue composition** (colorectal/endometrial/ovarian, unified by MMR-deficiency status rather than organ) is a disclosed deviation from the pilot's single-tissue design, per `report_GSE146889.md` Section 7; this may itself contribute additional non-tumor/normal structure that the confound-attribution framework (built around a single binary label) does not fully capture.
- **Quartile stratification is not label-balanced** (Section 4) — Q1 is 98% normal and Q4 is 100% tumor, because the classifier separates the classes well along this exact projection. The within-stratum control here is weaker evidence than a label-balanced stratification would give.
- **Bootstrap fixed-null approximation** (Section 6), as in the pilot, holds the reference null distribution fixed across resamples; a fully joint bootstrap would be a stronger version of this control.
- **Null draw counts (300 per condition)** are fewer than the pilot's 500/2000 for computational-budget reasons; p-values at the 0.0033 floor (1/301) indicate zero null draws exceeded the observed value, consistent with the reported z-scores, but a finer-resolution p-value would need more draws.
- **Residualization removes only the linear component of the confound** — as in the pilot, a quadratic/nonlinear residual difference between classes remains (evidenced by residualized CV AUC=0.206, itself significant at p=0.005 against a permutation null), and is not addressed by this control.
- **Loop-identity comparison (Section 5) is a single-draw, single-feature observation** — with n=45 (pre) vs n=22 (post) loop participants, this is suggestive but not a formal statistical test of "the loop changed"; it is reported descriptively.
- **Single dataset generalization:** this audit covers GSE146889 only. Two further replication datasets (CPTAC-CCRCC, TCGA-LUAD) remain to be audited under this same protocol before a cross-dataset confound-attribution conclusion can be drawn for the full replication program.

## 9. Artifacts

- `table1_within_class_decomposition_GSE146889.csv` — mixed/tumor-only/normal-only max-H1 vs. both nulls
- `table2_within_stratum_control_GSE146889.csv` — confound-quartile binning results
- `table3_residualization_control_GSE146889.csv` — pre/post residualization max-H1, z-scores, CV AUC, loop composition
- `table4_block_bootstrap_ci_GSE146889.csv` — bootstrap CI on mixed/tumor-only/normal-only signal-beyond-null
- `table5_confound_spectrum_GSE146889.csv` — PC-vs-confound correlation tiering
- `within_class_null_comparison_GSE146889.png` — Fig. 1 (mixed / tumor-only / normal-only vs Gaussian null)
- `within_stratum_and_residualization_GSE146889.png` — Fig. 2 (quartile z-scores; residualization z/AUC tradeoff)
- `bootstrap_ci_and_auc_check_GSE146889.png` — Fig. 3 (bootstrap CI; top-loop tumor-enrichment pre/post residualization)
