# Ablation / Hyperparameter-Sensitivity Sweep: PCA-Topology Cross-Domain Replication (GSE81089)

**Author:** Claude Science (TOPOLOGICA project) | **Date:** 2026-07-08 | **Dataset:** GEO [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089) (NSCLC bulk RNA-seq, refetched and reprocessed independently for this sweep)

## 1. Scope and framing (per universal-ablation-engine anti-pattern #8)

This sweep tests **robustness of the pre-registered finding** — "does PCA(50) inflate max-H1-persistence
beyond what a matched null model shows?" — to the pipeline's hyperparameter choices. Per the
universal-ablation-engine's own anti-pattern list (#8: "ablating hyperparameters as components — k=30
vs k=50 is a sweep, not ablation"), **items 1, 2, 4, and 5 below are hyperparameter sweeps, not
component ablations**: HVG gene count, PCA component count, and permutation count are continuous/ordinal
knobs, not discrete pipeline components with an on/off state. Only **item 3 (imputation rule)** is a
genuine Tier-1 component ablation — it compares discrete, qualitatively different rules for handling
missing data (zero-fill, mean-fill, median-fill, drop-affected-genes), which is the kind of
present/absent/alternative-implementation comparison Tier 1 classical ablation is built for.

Data were refetched from `ftp.ncbi.nlm.nih.gov` (`GSE81089_FPKM_cufflinks.tsv.gz`, 63,130 genes ×
218 samples) and reprocessed independently of any prior session state. The zero-fill center-point
configuration (HVG=2000, PCA=50, seed=42) reproduced the original pilot's real-data statistics exactly
(raw max-H1 = 3.887, PCA50 max-H1 = 4.671, PCA-delta = +0.784), confirming the refetch/reprocessing
pipeline is faithful to the pre-registered defaults before any sweep point is interpreted.

## 2. Deviations from the ideal protocol (compute-constrained, stated explicitly)

- **Permutation count reduced to n=500 for the gene-count, PC-count, imputation, and interaction sweeps**
  (vs. the pre-registered 2000). This was necessary for compute tractability: 18 unique configs × 2 null
  models × up to 4,000 dimensions each, on a 4-CPU sandbox, made 2000-permutation nulls per config
  prohibitive (single configs already took up to ~1 hour at n=500 with PC=100/HVG=4000). The
  permutation-count convergence check (Section 4) is exactly the safeguard against this deviation: it
  directly tests whether n=500 changes the conclusion relative to n=2000 at the pre-registered default,
  and finds it does not (z drifts by < 0.7 units of a ~15-24 unit statistic).
- **Gaussian null draws kept at the pre-registered 500** throughout — no deviation here.
- **Interaction check limited to the 2×2 corners** (HVG∈{500,4000} × PC∈{10,100}) rather than a full grid,
  per the universal-ablation-engine's own N≤5 full-factorial / corner-check guidance for a 2-factor design.
- **Mean/median imputation are slow** (~80-90s to impute the full 218×45,047 matrix via per-gene
  `fillna`, vs. <1s for zero-fill/drop) — absorbed into the per-config compute budget, not skipped.
- No config failed to compute; all 18 pre-planned sweep points completed and are reported (see Section 8
  for the full per-config compute-time ledger; nothing was silently dropped).

## 3. Complete sweep table

18 unique configurations (center point at HVG=2000/PC=50/zero-fill/n_perm=500 is shared across sweep
items 1, 2, 3, and 4, so it is computed once and reused — this is why 4 sweep dimensions produce 18,
not 4+4+4+5+4=21, unique runs).

| Config (tag) | Imputation | HVG | PCs | n_perm (pipeline) | Real raw max-H1 | Real PCA50 max-H1 | Real PCA-Δ | z (raw vs pipeline-null) | z (PCA vs pipeline-null) | z (raw vs Gaussian) | z (PCA vs Gaussian) | Pipeline-null Δ (mean) | Gaussian-null Δ (mean) | H1 holds (null-Δ ≥ real-Δ, both nulls) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HVG500_PC50 | zero | 500 | 50 | 500 | 3.265 | 2.096 | **−1.169** | 7.63 | 1.76 | 20.41 | 1.44 | 0.059 | 0.806 | ✅ |
| HVG1000_PC50 | zero | 1000 | 50 | 500 | 4.111 | 3.509 | **−0.602** | 13.22 | 4.91 | 27.60 | 5.07 | 0.722 | 1.268 | ✅ |
| **HVG2000_PC50 (default)** | zero | 2000 | 50 | 500 | 3.887 | 4.671 | **+0.784** | 15.80 | 5.04 | 24.29 | 5.17 | 1.518 | 1.996 | ✅ |
| HVG4000_PC50 | zero | 4000 | 50 | 500 | 4.342 | 5.865 | **+1.524** | 20.35 | 4.69 | 28.22 | 4.12 | 2.547 | 3.010 | ✅ |
| HVG2000_PC10 | zero | 2000 | 10 | 500 | 3.887 | 6.382 | **+2.495** | 15.80 | 13.41 | 24.29 | 14.03 | 1.092 | 1.444 | ❌ |
| HVG2000_PC25 | zero | 2000 | 25 | 500 | 3.887 | 5.163 | **+1.276** | 15.80 | 7.18 | 24.29 | 6.98 | 1.430 | 1.867 | ✅ |
| HVG2000_PC100 | zero | 2000 | 100 | 500 | 3.887 | 3.741 | **−0.146** | 15.80 | 2.86 | 24.29 | 3.02 | 1.354 | 1.768 | ✅ |
| HVG2000_PC50_mean | mean | 2000 | 50 | 500 | 4.325 | 4.914 | **+0.589** | 16.56 | 5.82 | 28.04 | 5.95 | 1.551 | 1.996 | ✅ |
| HVG2000_PC50_median | median | 2000 | 50 | 500 | 4.325 | 4.918 | **+0.593** | 18.73 | 6.09 | 28.04 | 5.96 | 1.518 | 1.996 | ✅ |
| HVG2000_PC50_drop | drop | 2000 | 50 | 500 | 4.302 | 4.734 | **+0.432** | 16.11 | 5.23 | 27.84 | 5.37 | 1.496 | 1.996 | ✅ |
| HVG2000_PC50_np100 | zero | 2000 | 50 | 100 | 3.887 | 4.671 | +0.784 | 15.71 | 5.25 | 24.29 | 5.17 | 1.483 | 1.996 | ✅ |
| HVG2000_PC50_np200 | zero | 2000 | 50 | 200 | 3.887 | 4.671 | +0.784 | 15.92 | 4.92 | 24.29 | 5.17 | 1.502 | 1.996 | ✅ |
| HVG2000_PC50_np1000 | zero | 2000 | 50 | 1000 | 3.887 | 4.671 | +0.784 | 15.52 | 5.14 | 24.29 | 5.17 | 1.520 | 1.996 | ✅ |
| HVG2000_PC50_np2000 | zero | 2000 | 50 | 2000 | 3.887 | 4.671 | +0.784 | 15.20 | 5.12 | 24.29 | 5.17 | 1.522 | 1.996 | ✅ |
| HVG500_PC10 (interaction corner) | zero | 500 | 10 | 500 | 3.265 | 4.627 | **+1.363** | 7.63 | 17.53 | 20.41 | 18.29 | −0.191 | 0.522 | ❌ |
| HVG500_PC100 (interaction corner) | zero | 500 | 100 | 500 | 3.265 | 2.548 | **−0.717** | 7.63 | 4.79 | 20.41 | 5.42 | −0.110 | 0.599 | ✅ |
| HVG4000_PC10 (interaction corner) | zero | 4000 | 10 | 500 | 4.342 | 9.994 | **+5.653** | 20.35 | 18.51 | 28.22 | 18.52 | 1.956 | 2.215 | ❌ |
| HVG4000_PC100 (interaction corner) | zero | 4000 | 100 | 500 | 4.342 | 6.787 | **+2.445** | 20.35 | 6.64 | 28.22 | 7.12 | 2.358 | 2.754 | ❌ |

Full table with p-values and standard deviations: [ablation_sweep_full_table.csv](ablation_sweep_full_table.csv).

**Threshold derivation (per Tier 1, no hardcoding):** sigma of the primary metric (real PCA-delta)
across all 18 configs = **1.469** (mean = 0.981). A hyperparameter effect is judged noise-floor below
0.5σ (0.735) and unambiguous above 2σ (2.939). By this data-derived yardstick, the swing in real
PCA-delta from −1.17 (HVG=500) to +5.65 (HVG=4000, PC=10) spans **4.64σ** — this is not a marginal or
noise-level sensitivity; hyperparameter choice materially changes the *sign and magnitude* of the
real-data PCA effect, even though (as shown below) the qualitative null-comparison conclusion mostly
survives it.

## 4. Permutation-count convergence: does n_perm=2000 suffice?

Bootstrap 95% CIs (2000 bootstrap resamples of the null distribution at each n_perm, per Tier 2)
on the z-score estimate itself, at the pre-registered default (HVG=2000, PC=50, zero-fill):

| n_perm | z (raw HVG), mean [95% CI] | z (PCA50), mean [95% CI] | Null mean (raw) | Null mean (PCA50) |
|---|---|---|---|---|
| 100 | 16.05 [12.72, 20.18] | 5.33 [4.47, 6.45] | 1.329 | 2.812 |
| 200 | 15.99 [13.88, 18.37] | 4.97 [4.38, 5.67] | 1.334 | 2.837 |
| 500 (used elsewhere in this sweep) | 15.85 [14.62, 17.25] | 5.05 [4.68, 5.45] | 1.318 | 2.836 |
| 1000 | 15.53 [14.70, 16.45] | 5.14 [4.89, 5.42] | 1.313 | 2.833 |
| 2000 (pre-registered) | 15.22 [14.63, 15.81] | 5.13 [4.94, 5.33] | 1.314 | 2.836 |

The null mean and standard deviation are essentially flat from n=100 through n=2000 (raw-space null
mean varies by <0.4% across the full range; PCA-space null mean by <0.9%). The z-score itself drifts by
only −0.63 (raw, 15.85→15.22) and +0.08 (PCA, 5.05→5.13) between the pre-registered n=500 and the
upper-bound n=2000 — both drifts are well within the bootstrap CI width at n=500, and the CI narrows
monotonically as n_perm increases (as expected from CLT scaling), rather than shifting location. Cohen's
d comparing the n=500 and n=2000 bootstrap z-distributions: d=−1.27 (raw), d=0.47 (PCA) — a
moderate-to-large *distributional* shift by Cohen's convention, but the practical drift in the point
estimate itself is under 1 unit against z-scores of 5–16.

**Verdict: n_perm=2000 (the pre-registered choice) has converged.** The permutation-count sweep gives no
evidence that going beyond 2000 would change the conclusion, and n=500 (used elsewhere in this sweep for
compute tractability) already gives an answer statistically indistinguishable from n=2000's.

Full table: [permutation_convergence_table.csv](permutation_convergence_table.csv).

## 5. Imputation-rule ablation (genuine Tier-1 component ablation)

Four imputation rules for the 1,124 sentinel (`-1.0`) values across 89 genes, all at the pre-registered
HVG=2000/PC=50/n_perm=500:

| Rule | Real raw max-H1 | Real PCA50 max-H1 | Real PCA-Δ | z (raw vs pipeline-null) | z (PCA vs pipeline-null) |
|---|---|---|---|---|---|
| zero-fill (pre-registered) | 3.887 | 4.671 | +0.784 | 15.80 | 5.04 |
| mean-fill | 4.325 | 4.914 | +0.589 | 16.56 | 5.82 |
| median-fill | 4.325 | 4.918 | +0.593 | 18.73 | 6.09 |
| drop-affected-genes | 4.302 | 4.734 | +0.432 | 16.11 | 5.23 |

**Output-set overlap (Jaccard), not just aggregate metric, per Tier-1 mandatory reporting:**

- **Top-H1-loop sample participation** (which of the 218 samples' vertices anchor the single most
  persistent H1 cocycle, extracted via `ripser`'s cocycle representatives): **Jaccard = 1.000 for every
  pair of imputation rules** — all four rules identify the *exact same 4 samples* (indices 14, 76, 132,
  170 in the sample ordering) as driving the top loop.
- **HVG-2000 gene-set overlap**: Jaccard ranges from 0.969 (zero-fill vs. drop) to 1.000 (mean vs.
  median), all ≥0.97 — imputation rule perturbs which genes enter the top-2000 HVG set by at most ~3%.

Full overlap table: [imputation_overlap_table.csv](imputation_overlap_table.csv). Imputation-rule summary:
[imputation_ablation_summary.csv](imputation_ablation_summary.csv).

**Verdict: imputation rule is inconsequential here.** All four rules keep z-scores well above the
3.0 replication threshold in both null comparisons (range 5.04–6.09 for PCA-space, 15.80–18.73 for
raw-space), all four give a positive, modest real PCA-delta (+0.43 to +0.78), and — most importantly —
all four converge on identical output structure (same 4-sample loop, near-identical HVG gene set). This
makes sense given the affected data is tiny relative to the full matrix (1,124 of 13.76M entries,
0.008%) — the ablation confirms this expectation rather than merely assuming it.

## 6. Interaction check: do gene-count and PC-count interact?

Per universal-ablation-engine Part 2 (joint removal vs. sum of individual LOO-style changes), using the
pre-registered center point (HVG=2000, PC=50, real PCA-delta=+0.784) as baseline:

| Corner | HVG-alone effect | PC-alone effect | Additive prediction | Observed joint effect | Interaction term | \|interaction\|/σ |
|---|---|---|---|---|---|---|
| HVG=500, PC=10 | −1.953 | +1.711 | −0.242 | +0.579 | **+0.820** | 0.56 |
| HVG=500, PC=100 | −1.953 | −0.930 | −2.882 | −1.501 | **+1.382** | 0.94 |
| HVG=4000, PC=10 | +0.740 | +1.711 | +2.451 | +4.868 | **+2.418** | 1.65 |
| HVG=4000, PC=100 | +0.740 | −0.930 | −0.190 | +1.661 | **+1.851** | 1.26 |

All four corners show an interaction term with magnitude between 0.5σ and 2σ of the sweep-derived
threshold — i.e., **ambiguous-to-moderate by the sigma-derived criterion, not clearly additive and not
clearly a strong independent interaction** at any single corner. However, the interaction term is
**positive in all four corners** (range +0.82 to +2.42), meaning the joint change is consistently larger
than the sum of the individual single-hyperparameter changes would predict — gene-count and PC-count
do not act independently on real-data max-H1-persistence; there is a systematic superadditive
component, most pronounced at the HVG=4000/PC=10 corner (interaction = 2.42, the only corner exceeding
1.5σ).

**Verdict: gene-count and PC-count interact, and the interaction is directionally consistent
(superadditive) rather than noise.** A one-hyperparameter-at-a-time sensitivity reading (Sections 3's
marginal sweeps alone) would understate how extreme the real-data statistic becomes at the true corner
of low-HVG/low-PC or high-HVG/low-PC space; the interaction check is not optional here per the
ablation-engine's own protocol, and this analysis found a genuine (if not overwhelming) interaction
requiring escalation flagging, not a clean pass.

Full table: [interaction_check_table.csv](interaction_check_table.csv).

## 7. Honest verdict: does the H1 finding hold robustly?

**Two distinct claims must be kept separate, per the pre-registration's own Section 1 framing:**

**H0 (does real data significantly exceed both nulls at all) — holds in 16/18 configs (89%).**
The two exceptions: HVG=500/PC=50 fails on **both** PCA-space null comparisons
(PCA-vs-Gaussian z=1.44, PCA-vs-pipeline z=1.76, both `<3.0`) — at only 500 genes, raw-space still passes
overwhelmingly (z=20.4) but PCA-space signal is weak against either null; and HVG=2000/PC=100 fails
marginally on PCA-space-vs-pipeline-null only (z=2.86, just under 3.0, while PCA-vs-Gaussian z=3.02
narrowly passes). Both failures occur in **PCA-reduced space specifically**, at parameter extremes (small
HVG or large PC-count relative to sample count), never in raw-feature space, where every single one of the
18 configs passes comfortably (z≥7.6, most ≥13).

**H1 (does PCA inflate real persistence LESS than or comparably to null inflation — the specific
pre-registered claim under test) — holds in 14/18 configs (78%), and fails in a hyperparameter-coherent
pattern, not scattered noise.** All 3 failures at PC=10 (regardless of HVG: 500, 2000, or 4000) show
real PCA-delta exceeding both null PCA-deltas — i.e., at very low PC-count, PCA inflation in the *real*
data outpaces the null-inflation artifact. The 4th failure (HVG=4000, PC=100) fails only the
pipeline-null comparison, not the Gaussian one. **PC=10 is a genuine boundary of the pre-registered
claim, not an artifact of this sweep's np=500 deviation** — it replicates across all three HVG levels
tested at that PC count. Its magnitude relative to σ differs by which analysis it is read from: as a
marginal deviation from the center point, HVG2000_PC10 alone is 1.16σ (delta_from_center/σ); as an
interaction-corner term (joint effect minus the additive prediction from single-hyperparameter changes),
HVG500_PC10 is 0.56σ and HVG4000_PC10 is 1.65σ — i.e. only the HVG4000_PC10 corner shows a large
interaction term in its own right, while HVG500_PC10's failure is driven mostly by the PC-alone effect
rather than an interaction with gene count.

**Overall: the H1 finding is directionally robust at the pre-registered PC=50 (all three tested gene
counts pass), degrades at the tested extremes, and fails in two distinct patterns rather than one.**
The dominant pattern is at low PC-count (PC=10): all three gene counts tested there (500, 2000, 4000)
show real PCA-delta exceeding both nulls, the opposite of the pre-registered directional claim, and this
replicates cleanly across gene count — the pre-registered default (PC=50) sits safely inside the region
where H1 holds (real PCA-delta +0.784, both nulls above it at +1.52 and +2.00), and the sweep's marginal
PC-count series at HVG=2000 (Section 3) shows the H1-supporting margin *shrinks monotonically* from
PC=25→50 and then reverses sharply at PC=10. This low-PC pattern is consistent with a known geometric
phenomenon: very-low-dimensional PCA embeddings are more prone to spurious loop formation from geometric
artifacts of the projection itself, independent of whether the underlying data is structured or noise —
precisely the OQ-006 "PCA fabricates phantom topology" mechanism the pilot analysis flagged, evidently
strongest at aggressive dimensionality reduction. **A second, distinct failure occurs at HVG=4000,
PC=100** — not a low-PC-count case — where the H1 criterion (the null's own PCA-driven delta should
match or exceed the real delta) holds against the Gaussian null (mean null delta +2.754 > real +2.445)
but fails against the pipeline-symmetric null (mean null delta +2.358 < real +2.445), even though the
real signal remains highly significant against both nulls in the ordinary sense
(z_pca_vs_pipeline=6.64, z_pca_vs_gauss=7.12 — this config fails H1 on the delta-comparison criterion,
not on a low z-score). This means H1's boundary is not fully described as "holds at PC≥25, fails only at
PC=10": at high gene count (4000), a failure also appears at PC=100 against the pipeline-symmetric null's
delta specifically. We report the sweep's actual 4-failure pattern rather than collapsing it to the
simpler PC=10-only story.

## 8. Anti-patterns avoided / deviations acknowledged

- **Anti-pattern #8 (ablating hyperparameters as components) avoided by explicit framing** (Section 1):
  items 1/2/4/5 are labeled sweeps throughout, not ablations; only item 3 (imputation) is called an
  ablation.
- **Anti-pattern #4 (ignoring output overlap) avoided**: Section 5 reports Jaccard overlap on both
  sample participation and gene sets for the one genuine component ablation, not just the aggregate
  z-score/delta.
- **Anti-pattern #7 (LOO without interaction check) avoided**: Section 6 explicitly runs the joint-vs-additive
  interaction check per Part 2's protocol, rather than treating the marginal HVG and PC sweeps as if they were independent.
  It found a real (if moderate) interaction and reports it rather than assuming additivity.
  The interaction magnitude was derived from the sweep's own sigma, not asserted.
- **Anti-pattern #1 (hardcoded thresholds) avoided**: the "does this hyperparameter matter" threshold
  (Section 3, 0.5σ/2σ) is derived from the sweep's own metric distribution, per Tier 1's explicit
  instruction, not chosen post-hoc to fit a narrative.
- **Deviation acknowledged**: n_perm=500 used instead of the pre-registered 2000 for 15 of 18 sweep
  configs, justified by compute tractability and validated post-hoc by the dedicated convergence check
  (Section 4), which found no meaningful difference between n=500 and n=2000 answers.
- **No config was fabricated or silently dropped.** All 18 pre-planned configurations computed
  successfully; per-config compute times ranged from 336.7s (n_perm=100) to 3576.5s (HVG=500/PC=100,
  the single slowest config, ~60 min) — see [ablation_sweep_full_table.csv](ablation_sweep_full_table.csv)
  for the complete compute-time ledger. Total sweep wall-clock time: ~5.4 hours on 4 CPU workers (19,487s summed across all 18 configs).

## 9. Limitations

- **Interaction check limited to 4 corners**, not a full grid over all HVG×PC combinations; a finer grid
  (e.g., HVG∈{500,1000,2000,4000} × PC∈{10,25,50,100}, 16 cells) would map the interaction surface more
  precisely than the 2×2 corner check used here.
- **n_perm=500 for 15/18 configs** is a compute-driven deviation from the pre-registered 2000, mitigated
  but not eliminated by the convergence check (which only directly validates the effect at the single
  center-point config, not at every sweep point).
- **Single real-data observation per config** (as in the original pilot/report.md) — no confidence
  interval can be placed on the real-data statistic itself, only on the null distributions it is compared
  against.
- **Imputation ablation's near-perfect overlap is partly a consequence of the sentinel values' tiny
  scale** (0.008% of entries) — this ablation would likely show more separation on a dataset with a
  larger missingness fraction, and this near-null result should not be over-generalized to datasets with
  substantially more missing data.
- **HVG=500/PC=100 was the single slowest config** (~60 min) due to `ripser`'s cost scaling with the
  ratio of PC count to sample count when PC count approaches the sample count (218); this bounded how
  finely the PC-count axis could be swept without a coarser or GPU-accelerated approach.

## 10. Artifacts

- `ablation_sweep_full_table.csv` — complete 18-config sweep table (all raw statistics, p-values, z-scores, compute times)
- `permutation_convergence_table.csv` — n_perm convergence check with bootstrap CIs
- `interaction_check_table.csv` — gene-count × PC-count interaction analysis
- `imputation_overlap_table.csv` / `imputation_ablation_summary.csv` — imputation ablation Jaccard overlaps and statistics
- `sweep_sensitivity.png` — 3-panel figure: gene-count sweep, PC-count sweep, permutation-count convergence
- `interaction_and_imputation_overlap.png` — interaction-check bar chart and imputation Jaccard-overlap heatmap
