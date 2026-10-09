# Usage

This file documents the commands present in the repository. This is a command-line research program, not a web service. The Python functions and generated tables are **internal interfaces**, not a versioned public application programming interface (API). Principal component analysis (PCA) and one-dimensional homology (H1) are the primary methods. Cohorts include Clinical Proteomic Tumor Analysis Consortium (CPTAC) clear-cell renal cell carcinoma (CCRCC) and The Cancer Genome Atlas lung adenocarcinoma (TCGA-LUAD).

## How do I obtain a first successful run?

Use a Bash terminal on Linux with Git, Python 3.11, and the Python `venv` module installed:

```bash
git clone https://github.com/smaniches/pca-topology-crossdomain-replication.git
cd pca-topology-crossdomain-replication
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
sha256sum -c checksums.sha256
python reproduce.py
test -s reproduce_output/persistence_diagrams_real_vs_gaussian.png
test -s reproduce_output/barcode_and_null_distributions.png
```

The last two commands exit unsuccessfully if the expected figure is missing or empty. The Python script also prints `REPRODUCE: SUCCESS` only if its default subprocesses and presence checks pass. This path loads committed `.pkl` checkpoints, so it tests dependency installation and figure production **without** rerunning the original stochastic study. Do not unpickle checkpoint files supplied by an untrusted party.

## What does `reproduce.py` accept?

Run these from the repository root with the virtual environment activated:

| Option | Behavior observed in `reproduce.py` | Important cost or condition |
| --- | --- | --- |
| No option | Regenerates two GSE81089 pilot figures from bundled checkpoint files | Fast installation check; does not recompute pilot statistics |
| `--full` | Additionally downloads GSE81089 and runs `code/03_statistics_and_results_table.py` | Fixed full pilot Gaussian/permutation draw counts; `--n-draws` **does not** reduce them |
| `--replication` | Runs GSE146889, CPTAC-CCRCC, TCGA-LUAD RNA sequencing and TCGA-LUAD methylation drivers | External data downloads may be required |
| `--ablation` | Runs the GSE81089 preprocessing, configuration generation, parameter sweep, aggregation, and default verification | Worker pool and cached sweep outputs |
| `--confound` | Runs the four historical cohort confound-audit drivers | Historical residualized AUC and null-comparison caveats apply |
| `--all` | Enables full pilot, replication, ablation, and confound phases | Launches all heavy phases; do not use as a first-run check |
| `--n-draws N` | Caps null draws in replication, sweep, and confound drivers; default 20 | A reduced count changes the Monte Carlo distribution; not the locked full-count result |

The default pilot figures run **before** any selected extra phase. The top-level script returns a nonzero process code when a subprocess fails or a checked output is missing. However, some child scripts have their own error-handling conventions; inspect printed output and saved tables rather than treating all logs as a common schema.

Examples:

```bash
python reproduce.py --replication --n-draws 20
python reproduce.py --ablation --n-draws 20
python reproduce.py --confound --n-draws 20
```

The first two examples may fetch data. These are reduced-draw smoke tests and are **not** a rerun at the original locked Monte Carlo counts.

## How do I run a single cohort instead of every cohort?

The standalone replication drivers take `--n-gauss` and `--n-perm`; they default to the original full-count values. Each supports `--skip-assert` for investigation, which suppresses its hard numeric assertions and should not be used to declare success. Run from the repository root:

```bash
python code/replication_GSE146889/01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
python code/replication_CPTAC_CCRCC/01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
python code/replication_TCGA_LUAD/01_rnaseq_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
python code/replication_TCGA_LUAD/02_methylation_fetch_preprocess_and_results_table.py --n-gauss 15 --n-perm 15
```

Only run the cohort you need. The corresponding `code/replication_*/README.md` identifies each original data source and assertion tolerance. These drivers write `*_results_table*.csv` into their respective script directories rather than `reproduce_output/`.

GSE146889 reads a public Gene Expression Omnibus (GEO) gene-count table and caches the downloaded file. CPTAC-CCRCC uses the committed Proteomic Data Commons (PDC) log2-ratio matrix and labels under `code/replication_CPTAC_CCRCC/data/`, or fetches from PDC GraphQL when the selected cache is absent. TCGA-LUAD drivers use the Genomic Data Commons (GDC) API and locally cached expression/methylation data where available. Both TCGA and CPTAC drivers accept `--data-dir` to select a cache location.

The methylation branch uses `PCA(35)` on its selected samples and measures the counted reduced-space H1 statistic with a spectral distance rather than Euclidean distance. Those values should not be compared across metrics without the stated distinction.

## How are the sweep and historical confound controls run independently?

To build and execute the pilot parameter sweep, use the same order as `reproduce.py --ablation`:

```bash
python code/ablation_sweep/00_download_and_preprocess.py
python code/ablation_sweep/01_make_configs.py
python code/ablation_sweep/02_run_sweep.py --quick 20 --workers 2
python code/ablation_sweep/03_aggregate_results.py
python code/ablation_sweep/verify_default.py
```

The preprocessing stage creates `code/ablation_sweep/data/imputation_variants.pkl`. The config builder creates `code/ablation_sweep/configs.json`. The sweep writes `code/ablation_sweep/results/sweep_results.json` and one `nulldist_<tag>.pkl` per completed config. Aggregation writes CSV reports there. The sweep supports `--tags` for a comma-separated subset of configuration names and `--workers` to adjust its multiprocessing pool. **Inspect all checkpoint entries for `FAILED: true` before trusting its completion message.** An existing failed tag is skipped on resume because the resume logic treats all saved tags as done.

The historical confound scripts reside in `code/confound_attribution_audit/`, one per cohort. Their `--quick N` caps Gaussian, permutation, and bootstrap draws. They also accept `--data-dir` and `--results-dir`; the GSE146889 driver additionally accepts `--local-null-dir`. They generate comma-separated values (CSV) summaries and strata tables, **not** a corrected cross-validated analysis. See [LIMITATIONS.md](LIMITATIONS.md) before interpreting residualization values.

## How do I reproduce or inspect the newer CPTAC sensitivity test?

The newer class-conditional analysis uses the committed CPTAC data and verifies its two Secure Hash Algorithm 256-bit (SHA-256) digests before running. To execute its synthetic unit checks and a bounded pilot:

```bash
python -m unittest discover -s audit/cptac_resid_null -p 'test_*.py' -v
python audit/cptac_resid_null/experiment.py --mode pilot --draws 24 --seed 20261009 --out reproduce_output/cptac_resid_null_pilot.json
```

The `experiment.py` command accepts:

| Argument | Contract in the code |
| --- | --- |
| `--mode pilot` | Requires 2–32 null draws; 24 by default |
| `--mode confirm` | Requires exactly 499 null draws; default 499 |
| `--draws N` | Overrides pilot count inside its allowed range; cannot change the confirm count |
| `--seed N` | Seeds NumPy's generator; default 20261009 |
| `--out PATH` | Required JSON output path; parent directories are created |

The JavaScript Object Notation (JSON) output contains observed max-H1, every null draw, mean, sample standard deviation, descriptive z-score, plus-one empirical tail probability, seed, runtime, software versions, and input hashes. `confirm` is the program's **execution-mode name**; scientifically the result is a *post-preregistration sensitivity analysis*, not a new preregistered confirmation. The archived complete result and its limitations are in `results/cptac_resid_null_sensitivity_20261009/`.

On GitHub, `.github/workflows/cptac-residualization-null.yml` runs unit checks and a 24-draw pilot on relevant pull requests. The 499-draw mode runs only after manual workflow dispatch with `mode=confirm`, on a standard hosted runner. The pipeline limits execution time but does not establish an unconditional cost or runtime guarantee for every environment.

## How do I check the archived independent GSE146889 results?

The following **read-only** verification checks the committed reconciliation tables and their recorded input digest without rerunning its archived multi-shard computation:

```bash
python audit/gse146889_clean_lineage/test_integrity_guards.py
(cd audit/gse146889_clean_lineage/results && sha256sum -c SHA256SUMS)
python audit/gse146889_clean_lineage/verify_reconciliation.py --repo-root .
```

The full clean-lineage workflow has separate `prepare.py`, `worker.py`, and `aggregate.py` scripts with explicit command-line arguments and sharded intermediate files. Its preserved `audit/gse146889_clean_lineage/EXECUTION_WORKFLOW.yml` is an **archive**, not an active workflow under `.github/workflows/`. Do not run it as part of an ordinary documentation or dependency check.

A separate full-count CPTAC re-execution is archived under `experiments/cptac_full_null_20261009/`, with its own run records. That is a replay of the original experiment, distinct from the conditional-null sensitivity analysis.

## How are outputs verified and how do I troubleshoot?

The tracked baseline list is `checksums.sha256`. From the root, `sha256sum -c checksums.sha256` checks listed paths; `audit/gse146889_clean_lineage/results/SHA256SUMS` checks its own audit outputs. The GitHub `CI` workflow performs those checks and runs `python reproduce.py` with Python 3.11. The presence of a file in Git history alone does not verify its scientific interpretation.

| Symptom | Inspect or change |
| --- | --- |
| `python3.11` or `venv` not found | Install Python 3.11 and its virtual-environment package; rerun installation in a new `.venv` |
| `ModuleNotFoundError` | Confirm the shell activated `.venv` and `python -m pip install -r requirements.txt` completed |
| Missing pilot checkpoint or `pickle.load` error | Confirm repository checkout and `sha256sum -c checksums.sha256` pass; do not source replacement pickle files from unknown parties |
| `reproduce_output/*.png` absent | Check the failed child script and its printed output; no-flag mode should create the directory and two figures |
| GEO/PDC/GDC HTTP, request, or download error | The optional phase requires a reachable upstream source; the default figure path does not |
| Cohort assertions fail | Check input cohort, cached file, preprocessing deviation, software versions, and metric; do not disable assertions to manufacture a pass |
| Different null p-value after `--n-draws` or `--quick` | Reduced draws have a different empirical tail resolution; compare with the corresponding full-count recorded study, not the smoke test |
| Sweep prints `SWEEP COMPLETE` with incomplete results | Inspect `sweep_results.json` for `FAILED: true` and missing tags; remove or repair the specific failed entry before resuming |
| SHA-256 mismatch after documentation edits | The tracked-file checksum registry must be updated deliberately; compare Git changes before accepting any new hash |

For module rationale read [ARCHITECTURE.md](ARCHITECTURE.md) and [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md); for scientific and operational limits read [LIMITATIONS.md](LIMITATIONS.md).
