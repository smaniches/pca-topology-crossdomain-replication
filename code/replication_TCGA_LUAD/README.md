# TCGA-LUAD replication cohort -- reproduction code

Reproduces the 3rd pre-registered confirmatory test in the TOPOLOGICA
PCA-topology program, across TWO omics layers from the same TCGA-LUAD
(lung adenocarcinoma) cases: RNA-seq (116 samples) and DNA methylation
(36 samples, the subset of cases with both layers available).

Full methodology and results: [`../../results/replication_TCGA_LUAD/tcga_luad_report.md`](../../results/replication_TCGA_LUAD/tcga_luad_report.md).

## Layout

- `_common.py` -- shared helper functions (HVG selection, standardization,
  intrinsic-dimension gate, persistent homology via ripser, the spectral/
  graph-Laplacian distance metric, null-model statistics). Not runnable
  directly; imported by both numbered scripts.
- `01_rnaseq_fetch_preprocess_and_results_table.py` -- RNA-seq layer.
- `02_methylation_fetch_preprocess_and_results_table.py` -- methylation layer.

## Running

```bash
# RNA-seq layer -- full pre-registered run (N_GAUSS=500, N_PERM=2000).
# ~20-40 minutes on a single modern core.
python3 01_rnaseq_fetch_preprocess_and_results_table.py

# Methylation layer -- full pre-registered run. Each null draw also does a
# graph-Laplacian eigendecomposition (spectral distance), so this is slower
# per-draw than the RNA-seq layer despite the smaller sample count (36
# samples): ~30-60 minutes on a single core for the full run.
python3 02_methylation_fetch_preprocess_and_results_table.py

# Fast smoke tests (a few minutes each) -- verify the pipeline runs
# end-to-end and reproduce the deterministic real-data numbers exactly.
python3 01_rnaseq_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
python3 02_methylation_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
```

Requires: numpy, pandas, scipy, scikit-learn, ripser, requests (only needed
if re-fetching from GDC; see below).

## What it reproduces

**RNA-seq layer:**

| Quantity | Expected value | Tolerance |
|---|---|---|
| Real max-H1 persistence, raw top-2000-HVG | 8.327 | 1e-3 (exact; independent of null draw count) |
| Real max-H1 persistence, PCA(50) | 7.433 | 1e-3 (exact; independent of null draw count) |
| 5-fold CV AUC (tumor vs normal, PCA50) | 0.998 | 0.02 |

**Methylation layer:**

| Quantity | Expected value | Tolerance |
|---|---|---|
| Real max-H1 persistence, PCA(35, Section-5-mandated SPECTRAL metric) | 0.0 EXACTLY | 1e-9 |

The methylation-layer zero is a reported **null-space collapse**, not a bug:
this cohort's PCA(35) space falls in the "transitional" intrinsic-dimension
regime, which the locked pre-registration protocol's Section 5 requires be
evaluated under a graph-Laplacian spectral distance metric rather than raw
Euclidean distance; under that metric, the observed H1 persistence is
exactly zero. The Euclidean-PCA max-H1 (4.09, printed by the script and
included in the "transparency only" columns of
`tcga_luad_methylation_results_table.csv`) is nonzero and NOT part of the
primary test family, per the Section-5 gate.

Both scripts print each expected number alongside a PASS/FAIL verdict and
assert them (unless `--skip-assert`). The real-data numbers above are
exact/assertable regardless of `--n-gauss`/`--n-perm`; only the
null-distribution z-scores/p-values in the results CSVs are noisier under a
reduced draw count.

## Data source

Both scripts look for a local cache under `./data/` first
(`rna_fpkm_matrix.pkl` + `rna_sample_meta.csv` for the RNA-seq layer;
`meth_beta_matrix.pkl` + `meth_sample_meta.csv` for the methylation layer)
and use it if present -- but **this cache is NOT bundled in this
repository** (the raw matrices are ~55MB and ~140MB respectively, too large
to commit as plain-text-diffable git objects for two omics layers of one
of three replication cohorts). On a fresh clone, both scripts re-fetch
directly from the public GDC API on first run (documented explicitly in
each script's docstring and `fetch_from_gdc()` function):
```
https://api.gdc.cancer.gov/files   (manifest query, filtered to TCGA-LUAD,
    Primary Tumor / Solid Tissue Normal cases with BOTH tumor and matched
    normal, restricted to the case overlap between RNA-seq and methylation)
https://api.gdc.cancer.gov/data/{file_id}   (per-file download)
```

## Disclosed deviations (all printed at runtime; see docstrings and
`tcga_luad_report.md` Sections 4/5/8 for full detail)

1. **RNA-seq bottom-variance selection:** 5556/60660 genes have exactly
   zero variance across all 116 samples, exceeding the pre-registered
   2000-gene bottom-variance selection size. Bottom-variance genes are
   selected among the variance>0 pool only.
2. **Methylation imputation/logit domain collision:** beta-values are
   bounded in [0,1]; the blanket "missing -> 0" rule is incompatible with
   `logit(beta)` diverging at beta=0. Probes NaN in ALL samples are dropped
   first; remaining NaNs are imputed to 0, then all values are clipped to
   `[1e-6, 1-1e-6]` before the logit transform.
3. **Methylation PCA cap:** PCA(50) is infeasible at n=36 samples (SVD cap);
   `n_components = min(50, n_samples-1) = 35` is used for methylation only.
4. **Methylation Section-5 metric gate:** the PCA(35) space falls in the
   "transitional" intrinsic-dimension regime, requiring the spectral
   (graph-Laplacian eigenmap) distance metric as the primary comparison
   metric instead of raw Euclidean distance (see "What it reproduces" above).
