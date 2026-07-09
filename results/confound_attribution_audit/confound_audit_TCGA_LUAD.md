# Confound-Attribution Audit: Is the TCGA-LUAD H1 Signal a Tumor/Normal Class-Mean-Shift Artifact?

**Project:** TOPOLOGICA PCA-topology cross-domain replication | **Dataset:** TCGA-LUAD (lung adenocarcinoma, RNA-seq layer only)
**Protocol:** `confound-attribution-audit` skill, identical methodology to the GSE81089 pilot audit
**Data reused** from the already-committed TCGA-LUAD confirmatory analysis (`rna_fpkm_matrix.pkl`, `rna_sample_meta.csv`, 116 samples: 58 Primary Tumor / 58 Solid Tissue Normal, STAR-Counts FPKM, GENCODE v36). The RNA-seq preprocessing pipeline (log1p(FPKM), top-2000-by-variance HVG, standardize, PCA(50), seed=42) was reproduced independently in this audit and matched the original report's real max-H1 values to full numerical precision (raw HVG: 8.326534 vs. reported 8.327; PCA(50): 7.433323 vs. reported 7.433), confirming the reproduction is faithful before any audit computation began. Methylation is excluded per task scope — that layer's PCA-space real max-H1 was already found to be exactly 0 (no signal exists to attribute).

## 1. Named confound

**Tumor/normal class-mean-shift magnitude** — identical operationalization to the GSE81089 pilot: the unit vector `u = (mean_tumor − mean_normal) / ‖mean_tumor − mean_normal‖` in standardized top-2000-HVG space, and each sample's scalar projection onto `u`. TCGA-LUAD's own confirmatory report already flagged this: CV AUC for tumor/normal prediction is 0.998 in both the raw HVG space and the PCA(50) space used for persistent homology — near-identical to GSE81089's 0.9995. This audit's job is the same as the pilot's: determine whether the detected H1 signal is this class-mean-shift in disguise.

## 2. Confound spectrum

| Quantity | Correlation with confound direction | Tier |
|---|---|---|
| PC1 | r = −0.9995, p = 6.5e-174 | **identical to the confound** |
| PC2–PC50 | \|r\| < 0.024, p > 0.80 (all) | free |
| 1-D confound-direction alone | — | AUC = 1.000 using only this single projection |
| Variance explained by the confound direction | 32.9% of total HVG variance (658.8 / 2000.0) | — |

**Reading:** the confound concentrates almost entirely in a single PC1 direction (r = −0.9995, effectively the confound itself), with PC2 onward completely free of it (all \|r\| < 0.024). This is an even cleaner (more extreme) separation than the pilot's, where the confound bled into PC2 and PC3 as well (r = −0.417 and −0.238). It also carries nearly 3× the variance share of the pilot's confound direction (32.9% vs. 11.2%) — the tumor/normal axis is a substantially larger fraction of total HVG variance in TCGA-LUAD than in GSE81089.

## 3. Frozen-baseline / group decomposition: within-class subsets vs. mixed set

| Condition | n | Observed max-H1 | z (Gaussian null) | p (Gaussian) | z (permutation null) | p (permutation) |
|---|---|---|---|---|---|---|
| Mixed (tumor+normal) | 116 | 8.327 | 58.73 | 0.0020 | 47.54 | 0.0005 |
| **Tumor-only** | 58 | **4.596** | **24.99** | 0.0020 | **20.60** | 0.0020 |
| Normal-only | 58 | 2.779 | 12.81 | 0.0020 | 11.18 | 0.0020 |

![Within-class null comparison]({{artifact:01934a5f-31d6-4fef-aca6-b0945cded736}})

**This does not replicate the pilot's decisive result in its exact form.** In GSE81089, the tumor-only subset reproduced the *identical* observed max-H1 value as the mixed set (3.887 = 3.887). Here, neither single-class subset comes close to reproducing the mixed-set magnitude: tumor-only (4.596) is 55% of the mixed value, and normal-only (2.779) is only 33%. Both single-class subsets clear their own null distributions by a wide margin (z=25.0 and z=12.8 against Gaussian, z=20.6 and z=11.2 against permutation — all far above the pre-registered z>3.0 bar), so **within-class topological structure independent of a mixed-class boundary genuinely exists in both tumor and normal tissue here**. But the *magnitude-matching* signature that made the pilot's finding so clean — "the loop's cocycle lives entirely inside one class and is exactly as large as the mixed-set loop" — is not present in TCGA-LUAD. The mixed-set max-H1 (8.327) exceeds both single-class values, meaning cross-class geometry (tumor-to-normal structure, or a genuinely larger topological feature only visible with both classes present) contributes materially to the mixed-set statistic in this dataset. This is a real, reportable difference in the confound-independence picture between the two datasets, not a reproduction failure of the audit itself (both nulls, both nested-space computations, and the raw-value cross-check all confirm the numbers).

## 4. Explicit control 1 — within-stratum effect (confound-quartile binning)

Tumor-only samples (n=58) were binned into quartiles of the confound projection (distance along the tumor/normal mean-shift direction), following the pilot's approach exactly.

| Quartile | n | Observed max-H1 | Null mean ± SD (Gaussian, matched n) | z | p (raw) | p (BH-adjusted) | Reject @0.05 |
|---|---|---|---|---|---|---|---|
| Q1 (least tumor-like) | 15 | 0.367 | 0.540 ± 0.235 | **−0.73** | 0.782 | 0.782 | **No** |
| Q2 | 15 | 1.524 | 0.547 ± 0.246 | 3.97 | 0.0040 | 0.0053 | Yes |
| Q3 | 14 | 2.676 | 0.515 ± 0.245 | 8.81 | 0.0020 | 0.0040 | Yes |
| Q4 (most tumor-like) | 14 | 1.359 | 0.495 ± 0.252 | 3.42 | 0.0020 | 0.0040 | Yes |

![Within-stratum and residualization]({{artifact:f4723762-6ab9-4eec-9c20-4dee19962eab}})

**This is a genuine partial non-replication of the pilot's within-stratum result.** In GSE81089 every quartile cleared the z>3.0 bar by a wide margin (min z=8.8). Here, **Q1 — the quartile of tumor samples least separated from normal along the confound axis — shows an observed max-H1 (0.367) that falls *below* its own null mean (0.540)**, giving a negative z-score and a BH-adjusted p=0.78, nowhere near significance. Q2–Q4 do clear the bar (z=3.4–8.8, BH-adjusted p<0.006 for all three), so the signal is not uniformly absent — it concentrates in the more "confound-typical" three-quarters of the tumor cohort and specifically fails in the quarter of samples that are least distinguishable from normal tissue along the confound direction. This is the opposite of what a pure-confound signal would predict (which would fail everywhere, or fail preferentially in the *most* confound-extreme quartile, not the least) but it is also not the clean "signal persists at every fixed confound value" result the pilot delivered. The honest reading: **the within-stratum control is mixed for TCGA-LUAD — 3 of 4 strata support confound-independence, 1 of 4 (n=15, the smallest and least tumor-like quartile) does not,** and quartile sample sizes here (14–15) are notably smaller than the pilot's (49–50), adding a real power caveat to Q1's null result that cannot be fully disentangled from a true absence of signal.

## 5. Explicit control 2 — residualization on the confound direction

The linear class-mean-shift direction was regressed out exactly (each sample's class-specific mean shifted to the pooled grand mean, per-gene, in standardized HVG space; class means matched to numerical precision post-residualization, max abs difference = 8.3e-16), then max-H1 persistence was recomputed in **both** the raw HVG(2000) space and the full HVG→PCA(50) pipeline space, against matched Gaussian and pipeline-symmetric permutation nulls.

| Space | Observed max-H1 | z (Gaussian) | z (permutation) | 5-fold CV AUC (tumor/normal) |
|---|---|---|---|---|
| Full HVG raw (class-mean intact) | 8.327 | 58.73 | 47.54 | 0.998 |
| **Class-mean-residualized HVG raw** | **2.779** | **14.48** | **12.13** | **0.000** |
| Full HVG→PCA(50) (class-mean intact) | 7.433 | 9.20 | 8.34 | 0.998 |
| **Class-mean-residualized HVG→PCA(50)** | **3.686** | **0.50** | **0.91** | **0.210** |

![Bootstrap CI and AUC check]({{artifact:fea53e4e-df9f-4571-b7e8-e3c56fc6e63f}})

**This is the most important departure from the pilot, and it is dataset-specific, not a generic audit failure.** Two results emerge depending on which space is tested:

- **Raw HVG (2000-d) space: the signal survives residualization convincingly** (z=14.48 Gaussian, z=12.13 permutation, both far above the z>3.0 bar), while the confound check confirms complete removal (CV AUC collapses from 0.998 to exactly 0.000 — below-chance, the same qualitative pattern as the pilot's residual-variance-structure finding, here even more extreme).
- **PCA(50) space — the space used throughout this dataset's own primary pre-registered test family and the pilot's headline comparison — the signal collapses to null-indistinguishable** (z=0.50 Gaussian, z=0.91 permutation, both far *below* the 3.0 bar) while CV AUC drops only partially (0.998 → 0.210, itself significantly different from chance under a label-permutation test, p=0.005, n=200 permutations — mirroring the pilot's own non-zero residual AUC finding, just at a different value).

A structural check clarifies why: because class-mean residualization is an exact per-class constant shift, within-class pairwise distances (and therefore within-class topology) are mathematically invariant under it — the residualized raw-space max-H1 (2.779) is **identical** to the pre-residualization normal-only value from Section 3 (2.779, verified to full floating-point precision). The dominant loop surviving residualization in raw space is the intra-normal-tissue loop, not a tumor-associated one — a different loop from the one carrying the mixed-set and tumor-only signal in Sections 3–4. Passing this control in raw space, then, does not straightforwardly confirm that the **tumor-associated** signal audited in Sections 3–4 is confound-independent; it confirms that *some* within-class topological structure (here, seemingly within the normal-tissue block) survives linear confound removal. In PCA(50) space — where the recombination of the 2000 residualized dimensions into 50 components evidently redistributes what little cross-class structure was needed to keep the loop's persistence above null — the signal is fully explained by the removed confound.

## 6. Interpreting the negative real-PCA-delta context

TCGA-LUAD's own confirmatory report found real PCA-delta = **−0.893** for RNA-seq (real max-H1 is *lower* in PCA(50) space than in raw HVG space: 7.433 < 8.327), while both matched nulls showed a strongly *positive* delta (+2.51 Gaussian, +2.35 permutation) — PCA inflates persistence far more on pure noise than on the real transcriptome. That earlier result answers a different question than this audit: it establishes that PCA does not *manufacture* the appearance of topological structure beyond what a matched null would show under the same dimensionality reduction. It says nothing about whether the *already-established* real signal is confound-independent — that is precisely this audit's job, and the two turn out to point in different directions here.

The negative real-PCA-delta means PCA(50) was already the *less generous* space for the real data even before any residualization — raw HVG carries more topological signal than its own PCA-reduced version. Combined with this audit's finding that PCA(50) space is exactly where the signal collapses under confound removal (Section 5), a coherent picture emerges: **whatever intra-tumor/intra-normal structure survives confound removal is concentrated in directions of the raw 2000-gene space that PCA(50) does not preserve.** PCA(50) here captures primarily (perhaps disproportionately) the tumor/normal class-mean-shift direction itself — consistent with the confound spectrum in Section 2 showing PC1 carries essentially the entire confound (r=−0.9995) and 32.9% of total variance, nearly 3× the pilot's confound-variance share. Once that dominant, confound-saturated direction is corrected for, the 50-dimensional PCA representation is left without enough independent topological content to clear null, even though the raw 2000-dimensional space still has it elsewhere. This is a stronger, more literal version of "PCA compresses away the confound-independent signal" than the pilot exhibited, and it means **the negative real-PCA-delta does not protect the confound-independence finding — if anything, it foreshadows exactly the PCA-space fragility this audit uncovers.**

## 7. Block-bootstrap CI

Bootstrap resampled the tumor block (n=58) with replacement, 500 resamples, recomputing observed max-H1 persistence each draw against the fixed tumor-matched Gaussian null (mean=0.867, SD=0.149, from the 500-draw reference null in Section 3).

| Statistic | Point estimate | 95% CI | Excludes 0? |
|---|---|---|---|
| Δ (observed max-H1 − Gaussian null mean) | 3.729 | [0.545, 4.989] | **Yes** |
| z (Δ / null SD) | 24.99 | [3.652, 33.430] | **Yes, but narrowly at the lower bound** |

The 95% CI excludes zero, but its lower bound (z=3.65) sits far closer to the pre-registered z>3.0 significance threshold than the pilot's did (pilot lower bound: z=11.0). This is a materially weaker, more fragile margin of safety than the pilot's bootstrap result — consistent with the smaller, more variable within-class signal already visible in Sections 3–4.

## 8. Verdict

Applying the skill's reading rubric to TCGA-LUAD's RNA-seq layer specifically, **the result is mixed, and materially weaker than the pilot — this audit does not reproduce the pilot's clean "genuine signal, fully confound-independent" finding.**

- **What survives:** Both the tumor-only and normal-only subsets show real topological structure clearing their own null distributions by wide margins (z=25.0 and z=12.8 against Gaussian null). Residualization in the **raw HVG (2000-d) space** leaves a signal intact (z=14.5 Gaussian / z=12.1 permutation) despite CV AUC collapsing to 0.000 — **but this raw-space residualized value (2.779) is numerically identical to the pre-residualization normal-only statistic from Section 3, confirming (per Section 5's structural check) it is the intra-normal-tissue loop passing through unchanged, not a demonstrated tumor-associated signal surviving the control.** The block-bootstrap CI on the within-tumor signal-beyond-null (a separate computation, Section 7) excludes zero on its own terms.
- **What does not survive, or replicates only partially:**
  - The tumor-only and normal-only subsets do **not** reproduce the mixed-set magnitude (unlike the pilot's exact match) — cross-class structure contributes materially to the mixed-set statistic here.
  - The within-stratum control **fails in one of four quartiles** (Q1, the least tumor-like tumor samples, n=15): observed max-H1 falls below its own null mean, BH-adjusted p=0.78. This is a genuine gap the pilot did not have (pilot: all four quartiles cleared z>8.8).
  - Most importantly: residualization in **PCA(50) space — the space this replication program's primary tests and the pilot's own headline comparison were run in — collapses the signal to null-indistinguishable** (z=0.50 Gaussian / z=0.91 permutation, both far below the z>3.0 bar), even though CV AUC only partially drops (0.998→0.210). The pilot's PCA(50)-space residualized signal, by contrast, survived (z=31.1 Gaussian / z=2.97 permutation).
  - The bootstrap CI's lower bound (z=3.65) sits much closer to the significance threshold than the pilot's (z=11.0), signaling a less robust margin.
- **Scope qualifier:** this audit is specific to TCGA-LUAD RNA-seq (n=116: 58 tumor / 58 normal). It does not address the methylation layer (already excluded per task scope, real PCA-space max-H1=0) or generalize beyond this dataset.

**Bottom line, stated plainly per the task's explicit instruction:** the GSE81089 pilot's confound-independence finding — that the H1 signal is genuine and not a proxy for the tumor/normal class-mean-shift — is a **dataset-dependent result, not a universal one**. In TCGA-LUAD, the signal shows real within-class structure and survives confound removal in raw feature space, but it does **not** cleanly survive the same audit in PCA(50) space (where the pre-registration's own tests are run), and it does **not** hold uniformly across every confound stratum. This is not a failure of the audit's execution — the underlying real-vs-null test (Section 7 of the original confirmatory report) still passes robustly in this dataset — but it is a genuine, honestly-reported non-replication of the specific "confound-independent in every tested sense" conclusion reached on the pilot dataset.

## 9. Limitations

- **Small quartile sizes.** Tumor-only quartiles here (n=14–15 each) are notably smaller than the pilot's (n=49–50), reducing power to detect a within-stratum signal and making Q1's null result harder to interpret cleanly as "absent" versus "underpowered." This caveat cuts against over-interpreting Q1's failure as definitive, but it does not remove the finding — Q1's observed statistic falling *below* its own null mean is not simply "non-significant," it is directionally opposite to what a genuine signal would show.
- **Residualization removes only the linear (mean-shift) component.** As in the pilot, a quadratic or higher-order component of tumor/normal separability is not addressed (evidenced by the non-chance residual AUC of 0.210 in PCA space and 0.000, i.e. below-chance, in raw space — both significant against a permutation null, p=0.005).
- **The raw-space residualized signal is not confirmed to be the same loop as the tumor-associated signal audited in Sections 3–4.** Because class-mean residualization is an exact constant per-class shift, within-class topology (and hence the raw-space residualized max-H1) is mathematically identical to the pre-residualization single-class value — here traced to the normal-only subset's own internal structure, not necessarily the tumor-associated loop the mixed-set and tumor-only analyses in Sections 3–4 concern. This audit did not perform a full cocycle trace to confirm loop identity across sections (as the pilot did with explicit sample-level tracing); this is a natural follow-up.
- **Single confound, linear only.** As in the pilot, this audit tests one named confound (tumor/normal class-mean-shift) via one method (linear residualization); other confounds (batch, purity, subtype) are out of scope.
- **BH-FDR family.** The 4-quartile family and 3-condition within-class family were each corrected independently (BH, alpha=0.05) rather than pooled into one larger family; a stricter combined correction would not change any conclusion here (Q1 remains non-significant either way; Q2–Q4 and the within-class results remain significant either way) but is noted for completeness.

## 10. Artifacts

- `table1_within_class_decomposition_TCGA_LUAD.csv` — mixed/tumor-only/normal-only max-H1 vs. both nulls
- `table2_within_stratum_control_TCGA_LUAD.csv` — confound-quartile binning results, with BH-FDR correction
- `table3_residualization_control_TCGA_LUAD.csv` — pre/post residualization max-H1, z-scores, CV AUC (raw and PCA50 spaces)
- `table4_block_bootstrap_ci_TCGA_LUAD.csv` — bootstrap CI on the within-tumor signal-beyond-null
- `table5_confound_spectrum_TCGA_LUAD.csv` — PC-vs-confound correlation tiering
- `within_class_null_comparison_TCGA_LUAD.png` — Fig. 1 (mixed / tumor-only / normal-only vs both nulls)
- `within_stratum_and_residualization_TCGA_LUAD.png` — Fig. 2 (quartile z-scores; residualization survival by space)
- `bootstrap_ci_and_auc_check_TCGA_LUAD.png` — Fig. 3 (bootstrap CI; AUC collapse confirmation)
- `tcga_luad_confound_audit_checkpoint.pkl` — preprocessed matrices, labels, confound projections, all null distributions, bootstrap draws (checkpoint; full null generation totaled ~10 minutes of compute)
