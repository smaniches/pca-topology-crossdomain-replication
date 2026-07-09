# Confound-Attribution Audit: Is the CPTAC-CCRCC H1 Signal a Tumor/Normal Class-Mean-Shift Artifact?

**Project:** TOPOLOGICA PCA-topology-cross-domain replication | **Dataset:** CPTAC-CCRCC proteomics (PDC000127)
**Protocol:** `confound-attribution-audit` skill, replicating the GSE81089 pilot's exact methodology on a second dataset (proteomics, renal cell carcinoma)
**Data reused** from the already-saved, pipeline-verified artifacts `cptac_ccrcc_log2ratio_raw.csv` (9,591 proteins × 194 samples, log2-ratio-to-reference, PDC's own median-centered TMT output) and `cptac_ccrcc_labels.csv` (110 Primary Tumor / 84 Solid Tissue Normal). Preprocessing exactly reproduced from `cptac_ccrcc_report.md`: impute missing → 0, median-center each sample (log-ratio space, confirmed a no-op — PDC's log2_ratio arrives pre-centered), top-2000 HVG by variance, per-protein standardization, PCA(50), seed=42. Reproduction verified: raw top-2000-HVG max-H1 persistence = 3.832947 and PCA(50) max-H1 persistence = 4.683274, both matching the report's stated values (3.833, 4.683) to 6 decimal places, and 5-fold CV AUC on PCA(50) = 1.000, matching the report's stated AUC exactly.

## 1. Named confound

**Tumor/normal class-mean-shift magnitude** — the direction and distance separating the two classes' means in standardized HVG(2000) space, i.e. the exact linear quantity a logistic-regression classifier exploits to achieve near-ceiling AUC. Operationalized identically to the GSE81089 audit: unit vector `u = (mean_tumor − mean_normal) / ‖mean_tumor − mean_normal‖` in standardized top-2000-HVG feature space, and each sample's scalar projection onto `u`.

CPTAC-CCRCC reports 5-fold CV AUC = 1.000 for tumor/normal prediction from PCA(50) — the same near-ceiling separability pattern seen in all four replication datasets in this program. This audit determines whether the detected H1 signal is this same class-mean-shift in disguise, or genuine structure independent of it.

## 2. Confound spectrum

| Quantity | Correlation with confound direction | Tier |
|---|---|---|
| PC1 | r = −0.9995, p ≈ 6.3e-290 | **identical to the confound** |
| PC2 | r = −0.027, p = 0.705 | free |
| PC3 | r = 0.014, p = 0.850 | free |
| PC4–PC50 | \|r\| < 0.005, p > 0.97 (all) | free |
| PC1 alone | — | AUC = 0.997 |
| 1-D confound-direction alone | — | AUC = 0.998 |
| Variance explained by the confound direction | 41.0% of total standardized HVG variance | — |

**Reading:** as in the GSE81089 pilot, the confound is concentrated almost entirely in a single dominant projection (here, PC1 itself — even more concentrated than the pilot, where the confound spanned PC1+PC2+PC3). PC2 onward are effectively free of it (|r| < 0.03). This dataset shows a *stronger* single-axis concentration of the confound than GSE81089, and the confound direction alone captures a larger share of total variance (41.0% vs. 11.2% in the pilot) — consistent with CPTAC's near-perfect (AUC=1.000, vs. the pilot's very-high-but-not-perfect) class separability.

## 3. Frozen-baseline / group decomposition: within-class subsets vs. mixed set

Max-H1 persistence computed separately in the **mixed** set (tumor+normal, n=194, matches the primary CPTAC replication result), the **tumor-only** subset (n=110), and the **normal-only** subset (n=84), each tested against its own size-matched null (300-draw Gaussian, 300-permutation pipeline-symmetric).

| Condition | n | Observed max-H1 | z (Gaussian null) | p (Gaussian) | z (permutation null) | p (permutation) |
|---|---|---|---|---|---|---|
| Mixed (tumor+normal) | 194 | 3.833 | 24.25 | 0.0033 | 21.57 | 0.0033 |
| **Tumor-only** | 110 | **3.833** | **18.60** | 0.0033 | **17.14** | 0.0033 |
| Normal-only | 84 | 2.319 | 9.20 | 0.0033 | 14.25 | 0.0033 |

![Within-class null comparison]({{artifact:6b3dd7ed-d9c2-4db8-abe6-c4f7b9eb3fc2}})

**This is the decisive first result, and it replicates the pilot's finding exactly.** The tumor-only subset reproduces the *exact same* observed max-H1 persistence value as the full mixed set (3.832947 in both, identical to 6 decimal places), with strong significance against both null models (z=18.6 Gaussian, z=17.1 permutation). Tracing the dominant loop's cocycle representative confirms this directly: the dominant H1 loop in the mixed set and in the tumor-only subset both have **the same three participating samples** (CPT0025880003, CPT0006440003, CPT0001540009), **all three of which are Primary Tumor** — the loop lives entirely inside one class, not on a tumor/normal boundary, exactly the pattern found in GSE81089's four-sample loop.

**Unlike the GSE81089 pilot** (where the normal-only arm, n=19, was too underpowered to draw any conclusion), CPTAC-CCRCC's normal-only subset (n=84) is well-powered and **also shows signal clearing the pre-registered z>3.0 bar** (z=9.2 Gaussian, z=14.2 permutation) — though at a lower observed max-H1 (2.319) than the tumor-only subset (3.833). This is a new result this audit can report that the pilot could not: in this dataset, normal tissue *also* carries topological structure beyond null, just weaker than the tumor-only signal. This does not change the audit's core question (whether the *observed effect* is the confound), but it means CPTAC-CCRCC lets us say something GSE81089 couldn't — the H1 signal is not exclusively a tumor phenomenon here, even though tumor carries the dominant loop.

## 4. Explicit control 1 — within-stratum effect (confound-quartile binning)

Tumor-only samples (n=110) were binned into quartiles of the confound projection itself (distance along the tumor/normal mean-shift direction), testing whether the H1 signal concentrates in one narrow band of "how tumor-like" a sample is, versus persisting at every fixed confound value.

| Quartile | n | Observed max-H1 | Null mean ± SD (Gaussian, matched n) | z | p |
|---|---|---|---|---|---|
| Q1 (least tumor-like) | 28 | 1.951 | 0.749 ± 0.202 | 5.97 | 0.0033 |
| Q2 | 27 | 1.681 | 0.745 ± 0.203 | 4.60 | 0.0033 |
| Q3 | 27 | 2.002 | 0.748 ± 0.203 | 6.17 | 0.0033 |
| Q4 (most tumor-like) | 28 | 3.033 | 0.749 ± 0.202 | 11.29 | 0.0033 |

![Within-stratum and residualization]({{artifact:bb2f0141-4215-42a0-a73b-02e4944626e5}})

**Every quartile clears the pre-registration's z>3.0 effect-size bar** (min z=4.60 at Q2, all four ≥4.60). This replicates the pilot's within-stratum result directionally (every quartile survives), though the margin here is narrower than in GSE81089 (which ranged z=8.8–28.4) — CPTAC-CCRCC's weakest quartile (Q2, z=4.60) is closer to the threshold than the pilot's weakest quartile (Q1, z=8.81). There is also a monotonic-looking increase toward Q4 (the most tumor-like quartile shows the largest z and largest observed statistic), a pattern not seen as clearly in the pilot — worth flagging as a partial confound-gradient rather than a fully flat stratification. Even so, **the signal does not collapse in any stratum**, including Q1 and Q2, the tumor samples least separated from normal along the confound axis — the core within-stratum-effect control result holds.

## 5. Explicit control 2 — residualization on the confound direction

The linear class-mean-shift direction was regressed out exactly (each sample's class-specific mean shifted to the pooled grand mean, per-protein, in standardized HVG space), then the identical PCA(50) pipeline was re-fit on the residualized matrix and max-H1 persistence recomputed against both a matched Gaussian null (dimension 50, 300 draws) and a pipeline-symmetric permutation null (raw-space permutation → re-fit PCA(50), 300 permutations).

| Space | Observed max-H1 | z (Gaussian) | z (permutation, pipeline-matched) | 5-fold CV AUC (tumor/normal) |
|---|---|---|---|---|
| Full HVG→PCA50 (class-mean intact) | 4.683 | 39.27 | 5.00 | 1.000 |
| **Class-mean-residualized HVG→PCA50** | **4.410** | **36.46** | **4.20** | **0.236** |

**Verification the residualization worked as intended:** class means are matched to numerical precision (max abs residual mean difference = 1.23e-15) — the dominant linear confound is verifiably gone. 5-fold CV AUC for tumor/normal collapses from 1.000 to 0.236 — below chance, echoing the GSE81089 pilot's own below-chance collapse (0.094), consistent with the residualized classifier's regularization picking up a residual noise structure rather than the removed linear signal.

**H1 signal survives, decisively.** Max-H1 persistence in the residualized PCA(50) space (4.410) barely moves from the intact space (4.683) and remains far outside both nulls (z=36.5 Gaussian, z=4.20 permutation) — both comfortably clearing the pre-registration's z>3.0 bar, and in fact *stronger* than the pilot's own post-residualization permutation-null margin (GSE81089: z=2.97, just under the z>3.0 bar; CPTAC-CCRCC: z=4.20, clearing it with room). **The topological structure is essentially unaffected by removing the exact quantity a linear classifier uses to separate the classes**, while that same classifier's performance is destroyed.

## 6. Block-bootstrap CI

Bootstrap resampled the **tumor block** (n=110, the natural group unit for this control) with replacement, 200 resamples, recomputing observed max-H1 persistence each draw and comparing against the fixed tumor-matched Gaussian null (mean=0.978, SD=0.154, held fixed across bootstrap draws per the pilot's efficiency-approximation convention).

![Bootstrap CI and AUC check]({{artifact:eb61cb5d-bbc4-4d6e-a38f-140ce897c81c}})

| Statistic | Point estimate | 95% CI | Excludes 0? |
|---|---|---|---|
| Δ (observed max-H1 − Gaussian null mean) | 2.855 | [0.833, 3.963] | **Yes** |
| z (Δ / null SD) | 18.60 | [5.43, 25.81] | **Yes** |

The 95% block-bootstrap CI on the within-tumor signal-beyond-null excludes zero by a wide margin — even the CI's weakest point (Δ=0.833, z≈5.4) clears the pre-registered z>3.0 bar. This is not a fragile point estimate, and the CI is comparably decisive to the pilot's own bootstrap result (GSE81089: Δ CI [1.37, 3.77], z CI [11.0, 30.4]).

## 7. Verdict

Applying the skill's reading rubric — **Δ CI excludes 0 and positive, survives the hard-negative/within-stratum control → genuine signal; credit it, name the specific feature/regime, do not inflate to "the method works."**

**GENUINE SIGNAL, independent of the tumor/normal confound, at the following precise scope — replicating the GSE81089 pilot's confound-independence finding in a second, independent dataset (different omics modality, different cancer type):**

- The H1 persistent-homology signal detected in CPTAC-CCRCC is **not** merely the tumor/normal class-mean-shift re-encoded as topology. The tumor-only subset alone reproduces the full mixed-set signal exactly (identical observed statistic to 6 decimal places, z=18.6 vs 24.2), the signal survives binning by the confound at every quartile (z=4.6–11.3, all clearing the pre-registered z>3.0 bar), and it survives exact linear removal of the confound direction (z=36.5 Gaussian / z=4.20 permutation, versus z=39.3/5.0 pre-residualization) while classifier AUC on the same residualized space collapses to 0.236 (below chance).
- Scope qualifier: this is demonstrated for the **tumor subpopulation** specifically (n=110) — the dominant H1 loop's three cocycle-participating samples are all tumor. **Unlike the pilot, this audit additionally shows the normal-only subset (n=84, well-powered here, unlike GSE81089's underpowered n=19) also clears the null threshold** (z=9.2 Gaussian, z=14.2 permutation) — a within-normal signal this program has not previously been able to test. This extends, rather than merely replicates, the pilot's finding: at least in CPTAC-CCRCC, confound-independent topological structure is not confined to the tumor class alone.
- The coexistence of near-ceiling CV AUC (1.000) alongside H1 signal beyond null in CPTAC-CCRCC is **not evidence of confounding** — the two are simultaneously true but statistically separable, exactly as found in GSE81089: the class-mean-shift explains the classifier's near-perfect AUC (unsurprising: the 1-D confound projection alone gives AUC=0.998), while a different, geometrically distinct part of the intra-tumor structure carries the topological signal.
- **This audit does not establish what the intra-tumor H1 loop biologically represents** — only that it is not a proxy for the tumor/normal label. A narrow tumor subtype, a batch effect, or genuine biological/proteomic heterogeneity within the tumor cohort all remain open explanations. Fisher's exact test on the loop's 3 tumor / 0 normal composition (reported in the original `cptac_ccrcc_report.md`, Section 9) gives p=0.260 — not statistically significant at this sample size, consistent with the underpowered-loop caveat noted there.

## 8. Proteomics-specific validity considerations

This is the **first application of the confound-attribution-audit protocol to a proteomics dataset** in this program (GSE81089 was bulk RNA-seq). Several proteomics-specific characteristics could plausibly affect the audit's validity and are flagged explicitly:

- **Missingness pattern.** 6.5% of entries in the raw log2-ratio matrix were non-detected/masked values, imputed to 0 per the pre-registration's rule. Proteomic missingness is frequently *not* missing-at-random — low-abundance proteins are systematically more likely to go undetected, and detection rates can differ between tumor and normal tissue for biological reasons (e.g., a protein overexpressed in tumor is more often *detected* in tumor samples and more often *missing→imputed-to-0* in normal samples). This means the confound direction itself could be partly an artifact of a tumor/normal difference in *missingness rate* rather than *quantitative abundance* — a distinct failure mode from anything the RNA-seq pilot faced (FPKM counts don't have the same structured non-detection problem at this scale). The residualization control (Section 5) removes the *linear mean* component of this effect regardless of its origin (biological or missingness-driven), so the survival of the H1 signal after residualization is robust to this concern for the *specific test performed* — but it means we cannot cleanly separate "genuine intra-tumor proteomic heterogeneity" from "residual structure induced by differential missingness patterns" as the origin of the surviving signal.
- **Dynamic range / log-ratio (not raw-intensity) input.** As disclosed in the original report, PDC's public API exposes only log2(ratio-to-common-reference), not linear intensities, requiring a deviation from the pre-registration's literal normalization instructions (median-centering in log-ratio space, skipping a second log2). This audit reused that already-validated deviation without modification; it does not introduce a new validity concern beyond what the original report already disclosed, but it means CPTAC-CCRCC's feature values have different scale/shape properties (signed, ratio-based) than GSE81089's (log1p of FPKM, non-negative) — a difference in what "distance" means geometrically between the two datasets that should temper any claim that the *magnitude* of signal (z-scores, persistence values) is directly comparable across datasets, even though the qualitative verdict (genuine, confound-independent) replicates.
- **Smaller dominant loop (3 vs. 4 participants).** The GSE81089 loop had 4 cocycle-participating tumor samples; CPTAC-CCRCC's has 3. Both are too few for a well-powered Fisher's exact test (both p>0.13/0.26), so this audit's within-class/within-stratum/residualization controls — which operate on the *aggregate* max-H1 statistic across the whole sample, not just the loop's own membership — remain the load-bearing evidence for "not the confound," not the loop-composition check, in both datasets.

## 9. Limitations

- **Fixed-null approximation in the bootstrap** (Section 6), as in the pilot, holds the reference null distribution fixed across resamples rather than regenerating it per-draw; justified by the null's much smaller SD (0.15) relative to the bootstrap spread of the observed statistic, but a fully joint bootstrap would be a stronger version of this control if compute allows.
- **Null draw count reduced from the pilot's 500/2000 to 300/300** (Gaussian/permutation) for compute-time reasons in this replication; z-scores at this draw count remain stable and well above threshold in every comparison, but exact p-values are capped at the resolution these draw counts allow (minimum achievable p ≈ 1/301 = 0.0033, versus the pilot's 1/501 or 1/2001).
- **Within-stratum quartile z-scores show a partial confound-aligned gradient** (Q4 > Q1,Q3 > Q2) not seen as clearly in the pilot's flatter stratification — while every quartile clears the threshold, this pattern is worth tracking in future datasets as a possible sign that *some* residual confound-correlated component contributes to signal strength even though it does not create the signal from nothing.
- **Loop-participant enrichment (Section 7 aside) is observational, not a formal test** — n=3 participating samples is too small to test statistically (Fisher's exact p=0.260 from the original report, not re-derived here).
- **Two datasets now audited** (GSE81089, CPTAC-CCRCC); whether this confound-independence result holds in the remaining two replication datasets (GSE146889, TCGA-LUAD) — both of which also showed the CV-AUC/H1-signal coexistence — remains untested and is the natural next step for this program.

## 10. Artifacts

- `table1_within_class_decomposition_CPTAC_CCRCC.csv` — mixed/tumor-only/normal-only max-H1 vs. both nulls
- `table2_within_stratum_control_CPTAC_CCRCC.csv` — confound-quartile binning results
- `table3_residualization_control_CPTAC_CCRCC.csv` — pre/post residualization max-H1, z-scores, CV AUC
- `table4_block_bootstrap_ci_CPTAC_CCRCC.csv` — bootstrap CI on the within-tumor signal-beyond-null
- `table5_confound_spectrum_CPTAC_CCRCC.csv` — PC-vs-confound correlation tiering
- `within_class_null_comparison_CPTAC_CCRCC.png` — Fig. 1 (mixed / tumor-only / normal-only vs both nulls)
- `within_stratum_and_residualization_CPTAC_CCRCC.png` — Fig. 2 (quartile z-scores; residualization survival)
- `bootstrap_ci_and_auc_check_CPTAC_CCRCC.png` — Fig. 3 (bootstrap CI; AUC collapse confirmation)
