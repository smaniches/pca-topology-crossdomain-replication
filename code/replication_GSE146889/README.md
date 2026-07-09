# GSE146889 replication cohort -- reproduction code

Reproduces the 2nd pre-registered confirmatory test in the TOPOLOGICA
PCA-topology program: does PCA(50) inflate persistent-homology (H1) signal
beyond a matched null, in an independent bulk RNA-seq cohort (colorectal /
endometrial / ovarian mismatch-repair-deficiency study, GSE146889, 176
samples)?

Full methodology and results: [`../../results/replication_GSE146889/report_GSE146889.md`](../../results/replication_GSE146889/report_GSE146889.md).

## Running

```bash
# Full pre-registered run (N_GAUSS=500, N_PERM=2000). Re-fetches the raw
# GEO count matrix on first run (~50MB download), then runs 2500 total
# persistent-homology computations. On the order of 30-60 minutes on a
# single modern core.
python3 01_fetch_preprocess_and_results_table.py

# Fast smoke test (a few minutes) -- verifies the pipeline runs end-to-end
# and reproduces the deterministic real-data numbers exactly, but the null
# distributions are too small to be a confirmatory result.
python3 01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
```

Requires: numpy, pandas, scipy, scikit-learn, ripser (see `../../requirements.txt`).

## What it reproduces

| Quantity | Expected value | Tolerance |
|---|---|---|
| Real max-H1 persistence, raw HVG (top-2000 by variance) | 6.204 | 1e-3 (exact; independent of null draw count) |
| Real max-H1 persistence, PCA(50) | 7.564 | 1e-3 (exact; independent of null draw count) |
| 5-fold CV AUC (tumor vs normal, PCA50 features) | 0.925 | 0.02 |
| Intrinsic dimension (MLE estimator, k=10), PCA(50) space | 8.53 | 0.1 (exact; independent of null draw count) |

The script prints each of these alongside the expected value and a PASS/FAIL
verdict, and asserts them (raising `AssertionError` on mismatch) unless
`--skip-assert` is passed. The real-data numbers above do not depend on
`--n-gauss`/`--n-perm` (they're computed from the fixed GEO download with
fixed seeds) and are asserted at full tolerance even in a reduced-draw smoke
test; only the null-distribution z-scores/p-values in
`final_results_table_GSE146889.csv` are noisier under a reduced draw count.

## Data source

Raw data is fetched directly from NCBI GEO's public FTP mirror on first run
(cached locally afterward as `GSE146889_GeneCount.tsv.gz` in this directory):

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/suppl/GSE146889_GeneCount.tsv.gz
```

No bundled/cached copy of the raw data is included in this repository (only
the derived `final_results_table_GSE146889.csv` output and figures under
`../../results/replication_GSE146889/` are committed) -- the fetch is fast
and the source is a stable public archive.

## Verified deviations from the pilot's literal preprocessing

None beyond what's already disclosed in the report -- see
`report_GSE146889.md` Section 2-4 for the one gene-row NaN-drop.
