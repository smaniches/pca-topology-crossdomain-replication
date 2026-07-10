# Ablation sweep (18 configurations, GSE81089 pilot)

Reproduces `results/ablation_sweep/ablation_sweep_full_table.csv` and the
other 4 supporting tables reported in `results/ablation_sweep/ablation_report.md`.

## What this tests

The core TOPOLOGICA finding is "does PCA inflate max-H1-persistence beyond
what a matched null model shows?" This sweep tests robustness of that
finding to the pipeline's hyperparameter choices on the GSE81089 NSCLC
pilot cohort (218 samples): HVG gene count, PC count, imputation rule for
the ~1124 `-1.0` sentinel ("missing") values in 89 genes, and permutation
draw count. Per the `universal-ablation-engine` skill's anti-pattern #8,
items 1/2/4/5 below are hyperparameter **sweeps**, not component ablations;
only item 3 (imputation rule) is a genuine discrete-component ablation.

18 unique configs result from the union of:
1. Gene-count sweep: HVG in {500, 1000, 2000, 4000} (PC=50, zero-fill, n_perm=500)
2. PC-count sweep: PC in {10, 25, 50, 100} (HVG=2000, zero-fill, n_perm=500)
3. Imputation-rule ablation: {zero, mean, median, drop} (HVG=2000, PC=50, n_perm=500)
4. Permutation-count convergence: n_perm in {100, 200, 500, 1000, 2000} (HVG=2000, PC=50, zero)
5. Interaction corners: HVG in {500, 4000} x PC in {10, 100} (zero-fill, n_perm=500)

The pre-registered default (HVG2000_PC50_zero_np500) is the shared center
point and is de-duplicated across items, giving 18 (not 21) configs.

## Files

- `00_download_and_preprocess.py` -- fetches GSE81089 FPKM data from GEO
  (`GSE81089_FPKM_cufflinks.tsv.gz`) and builds the 4 imputation-rule
  variants of the log1p/HVG-filtered expression matrix. Regenerates from
  scratch; does **not** read the repo's `results/pilot_GSE81089/checkpoints/`
  pickles (those are the original pilot run's artifacts, kept for reference).
- `01_make_configs.py` -- enumerates the 18-config grid to `configs.json`.
- `02_run_sweep.py` -- runs all 18 configs: for each, computes the real-data
  H1 persistence statistic in raw-gene space and PCA-projected space, and
  compares against both a pipeline-symmetric permutation null and a
  dimension-matched Gaussian null. Checkpoints after every config
  (`results/sweep_results.json`, resumable) and saves each config's full
  null distributions (`results/nulldist_<tag>.pkl`).
- `03_aggregate_results.py` -- builds the 5 final CSVs from
  `sweep_results.json` and the null-distribution pickles: the full sweep
  table, the permutation-convergence table (with bootstrap z-score CIs and
  Cohen's d for the n_perm=500-vs-2000 comparison), the interaction-check
  table (additive-vs-observed effect decomposition at the 4 HVG x PC
  corners), the imputation-overlap table (Jaccard similarity of the
  top-H1-loop cocycle vertex sets and HVG gene sets across imputation
  rules), and the imputation-ablation summary.
- `verify_default.py` -- checks the pre-registered default config's
  real-data statistics against the values committed in the repo's
  `ablation_sweep_full_table.csv`, to full float precision (these are
  deterministic given the data/PCA seed and do not depend on permutation
  count, so they hold exactly regardless of `--quick`).

## Data provenance

`00_download_and_preprocess.py` always re-fetches and reprocesses GSE81089
from the public GEO accession rather than reading the repo's pilot
checkpoints, so this sweep is independently reproducible from a fresh
clone with no prior state.

## Running

```bash
python 00_download_and_preprocess.py
python 01_make_configs.py
python 02_run_sweep.py --workers 4          # full run: ~hours (up to n_perm=2000)
# OR, for a quick correctness smoke-test (real-data stats unaffected, only
# null-distribution-derived z/p-values are approximate):
python 02_run_sweep.py --workers 4 --quick 20
python 03_aggregate_results.py
python verify_default.py
```

## Verification performed

Ran the full 18-config grid locally with `--quick 20` (20 permutation/
Gaussian draws per config instead of the pre-registered 500/2000, for
compute tractability on a 4-CPU sandbox -- this is a reduction from the
committed run and is disclosed here explicitly). Compared every config's
`real_raw_max_pers`, `real_pca_max_pers`, and `real_pca_delta` (all
deterministic given the data, preprocessing, and PCA seed -- independent of
permutation count) against the repo's committed
`results/ablation_sweep/ablation_sweep_full_table.csv`:

**All 18 configs matched to exact bit-for-bit precision (max absolute
difference = 0.0).** This confirms the extracted pipeline reproduces the
committed results exactly on the parts of the computation that do not
depend on the (stochastic, permutation-count-dependent) null distributions.
The pre-registered default config's `real_pca_delta = 0.7840385437011719`
and `z_pca_vs_pipeline` (which does depend on n_perm, so only checked
loosely under `--quick`) were verified via `verify_default.py`.

The 5 output tables' column structure (names and order) was verified to
match the committed CSVs exactly; row counts also matched (18 sweep rows, 5
convergence rows, 4 interaction rows, 6 overlap-pair rows, 4 imputation-
summary rows).

A full-precision run (`--quick` omitted, n_perm as specified per config)
was NOT executed in this session due to compute-time constraints (would
take multiple hours) -- this is the disclosed reduction. The deterministic
real-data statistics verified above do not depend on this choice.

## Figure generation

`04_figure_sweep_sensitivity.py` regenerates `paper/figures/fig3_ablation_sensitivity.png`
from the 4 committed CSVs above (`ablation_sweep_full_table.csv`,
`permutation_convergence_table.csv`, `interaction_check_table.csv`,
`imputation_overlap_table.csv`) -- no raw pickles or JSON needed. Recovered
via `host.lineage` from the original session (previously this figure existed
only as a static PNG with no committed generation script). Every plotted
value was checked against its source CSV before this script was written.
The regenerated PNG is NOT byte-identical to the committed one (SHA-256
differs); a pixel-level diff confirms the difference is confined to
0.12% of pixels (mean per-pixel channel delta 0.028 on a 0-255 scale;
isolated pixels reach a max delta of 111, consistent with anti-aliased
text/line edges shifting by a fraction of a pixel between matplotlib
renders, not a systematic color or data difference), consistent with
font-cache/rendering non-determinism across environments -- but this was
not confirmed against the exact rendering environment that produced the
originally-committed PNG.
