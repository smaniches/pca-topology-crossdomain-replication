# CPTAC-CCRCC replication cohort -- reproduction code

Reproduces the 1st pre-registered confirmatory test in the TOPOLOGICA
PCA-topology program, on an independent omics modality: does PCA(50)
inflate persistent-homology (H1) signal in renal cell carcinoma
proteomics (CPTAC-CCRCC, PDC study PDC000127, 194 samples)?

Full methodology and results: [`../../results/replication_CPTAC_CCRCC/cptac_ccrcc_report.md`](../../results/replication_CPTAC_CCRCC/cptac_ccrcc_report.md).

## Running

```bash
# Full pre-registered run (N_GAUSS=500, N_PERM=2000). Uses the cached
# preprocessed matrix under ./data/ (see "Data source" below); persistent
# homology uses gudhi's Vietoris-Rips complex. ~20-40 minutes on a single
# modern core.
python3 01_fetch_preprocess_and_results_table.py

# Fast smoke test (a few minutes) -- verifies the pipeline runs end-to-end
# and reproduces the deterministic real-data numbers exactly.
python3 01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
```

Requires: numpy, pandas, scipy, scikit-learn, gudhi, requests (only needed
if re-fetching from PDC; see below) (see `../../requirements.txt` --
`requests` is a transitive dependency of `GEOparse`, already listed there).

## What it reproduces

| Quantity | Expected value | Tolerance |
|---|---|---|
| Real max-H1 persistence, raw top-2000-HVG | 3.833 | 1e-3 (exact; independent of null draw count) |
| Real max-H1 persistence, PCA(50) | 4.683 | 1e-3 (exact; independent of null draw count) |
| 5-fold CV AUC (tumor vs normal, PCA50) | 1.000 | 0.02 |

The script prints each alongside the expected value and a PASS/FAIL verdict,
and asserts them unless `--skip-assert` is passed. The real-data numbers
above are exact/assertable regardless of `--n-gauss`/`--n-perm`; only the
null-distribution z-scores/p-values in `cptac_final_results_table.csv` are
noisier under a reduced draw count.

## Data source

By default, this script reads a cached, already-preprocessed copy of the
log2-ratio matrix from `./data/cptac_ccrcc_log2ratio_raw.csv` and
`./data/cptac_ccrcc_labels.csv` (bundled in this repository -- see
`checksums.sha256`). If that cache is absent, it re-fetches directly from
the public Proteomic Data Commons (PDC) GraphQL API:

```
https://pdc.cancer.gov/graphql
  query: quantDataMatrix(pdc_study_id: "PDC000127", data_type: "log2_ratio", acceptDUA: true)
  query: biospecimenPerStudy(pdc_study_id: "PDC000127", acceptDUA: true)
```

**Disclosed deviation:** the PDC public API only exposes CCRCC protein
abundance as log2(ratio-to-common-reference), not raw linear-scale
intensities. The pre-registration's "median-normalize across samples, then
log2 transform" rule (written assuming raw linear intensity input) is
therefore applied as median-CENTERING in log-space, skipping the redundant
log2 step. This is printed explicitly at runtime and documented in the
script's module docstring and in `cptac_ccrcc_report.md`.
