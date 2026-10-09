# Architecture

## What runs when a maintainer invokes the project?

`reproduce.py` is the top-level Python command-line orchestrator. It calls other Python scripts with `subprocess.run` and the active interpreter (`sys.executable`), checks exit codes, and checks for expected pilot output files. The cohort labels refer to the Gene Expression Omnibus (GEO), Proteomic Data Commons (PDC), Genomic Data Commons (GDC), Clinical Proteomic Tumor Analysis Consortium clear-cell renal cell carcinoma (CPTAC-CCRCC), and The Cancer Genome Atlas lung adenocarcinoma (TCGA-LUAD). RNA-seq means RNA sequencing, H1 means one-dimensional homology, and API means application programming interface.

The no-flag path creates `reproduce_output/` and runs `code/04a_figure_persistence_diagrams.py` and `code/04b_figure_barcode_null_distributions.py`. Both load committed pilot checkpoint pickles. This path redraws figures; it neither downloads new samples nor recomputes the original full pilot null distributions.

Optional flags launch other phases from `reproduce.py`. These phases have distinct data-fetching and output conventions, so there is no single normalized dataset loader or result schema.

```text
reproduce.py
  default --> pilot checkpoint pickles --> 04a/04b figure scripts
                 |                            |
                 +----------------------------+--> reproduce_output/*.png
  --full --> code/03_statistics_and_results_table.py
              --> GEO GSE81089 --> log1p/feature selection --> H1 and nulls
              --> reproduce_output/final_results_table.csv
  --replication --> code/replication_GSE146889/01_*.py --> GEO
                --> code/replication_CPTAC_CCRCC/01_*.py --> cached PDC or PDC API
                --> code/replication_TCGA_LUAD/01_*.py --> GDC RNA-seq
                --> code/replication_TCGA_LUAD/02_*.py --> GDC methylation
  --ablation --> code/ablation_sweep/00_*.py --> 01_*.py
             --> 02_*.py --> 03_*.py --> verify_default.py
  --confound --> code/confound_attribution_audit/confound_audit_{cohort}.py
                  --> confound_audit_common.py
```

The diagram describes dispatch, not shared in-memory execution. Each script is a new process and may reread data. The top-level driver prints `REPRODUCE: SUCCESS` only when its subprocess exit-code and selected file-presence checks pass; it is not an automatic validation of scientific assumptions.

## What are the analysis components?

**Pilot and figures.** `code/03_statistics_and_results_table.py` fetches the Gene Expression Omnibus (GEO) GSE81089 expression table, handles sentinel values, uses a `log1p` transform, selects genes by variance, computes principal component analysis (PCA) and Vietoris–Rips persistent homology, and writes a comma-separated values (CSV) table. The 04a/04b scripts instead read `results/pilot_GSE81089/checkpoints/*.pkl` and generate the image outputs. Their own Gaussian examples are display data; they do not replace the stored full null runs.

**Replications.** Separate drivers under `code/replication_GSE146889/`, `code/replication_CPTAC_CCRCC/`, and `code/replication_TCGA_LUAD/` own cohort download, data cleaning, null generation, and expected-number assertions. The Cancer Genome Atlas lung adenocarcinoma (TCGA-LUAD) branch has a shared `_common.py` with intrinsic-dimension diagnostics, persistence summaries, spectral-distance construction, and Benjamini–Hochberg multiple-testing utilities. The Clinical Proteomic Tumor Analysis Consortium clear-cell renal cell carcinoma (CPTAC-CCRCC) driver uses GUDHI; other core branches use ripser. Data and metric differences are not hidden behind a single interchangeable engine.

**Parameter study.** `code/ablation_sweep/00_download_and_preprocess.py` writes a pickle containing imputation variants; `01_make_configs.py` writes a JavaScript Object Notation (JSON) configuration array; `02_run_sweep.py` evaluates each configuration with a multiprocessing pool and checkpoints results after every configuration; `03_aggregate_results.py` writes derived CSV tables; `verify_default.py` checks the declared baseline. The sweep resumes by configuration tag and does not automatically retry tags previously written with an error record.

**Historical confound audit.** `confound_audit_common.py` provides within-class comparisons, confound-score quartiles, label-conditioned class-mean residualization, and row-bootstrap calculations. Cohort-specific drivers provide their own preprocessing. The residualized cross-validation area is a *known methodological limitation*; sharing the function does not establish predictive validity. See [LIMITATIONS.md](LIMITATIONS.md).

**Independent audit paths.** `audit/gse146889_clean_lineage/` has a data-preparation step, sharded Monte Carlo worker, aggregator, and read-only reconciliation verifier. It is separate from the transcript-reconstructed historical GSE146889 driver. `audit/cptac_resid_null/experiment.py` verifies bundled CPTAC files by Secure Hash Algorithm 256-bit (SHA-256), shuffles each protein within each label class, applies the same class-mean residualization and PCA to observed and surrogate inputs, then computes max-H1 with ripser. Its `test_experiment.py` contains synthetic and numerical checks. The standalone `experiments/cptac_full_null_20261009/` records a full-count replay of the *original* CPTAC analysis, not a new cohort. A second distinct `audit/cptac_covariance_null/` experiment takes tumor-only samples from the same cached preprocessed CPTAC matrix and uses orthogonal rotations of the sample-centered row subspace. Its exact singular-value decomposition retains the complete feature-covariance spectrum, and its recorded 499-draw result lives under `results/cptac_tumor_covariance_null_20261009/`. It tests a different conditional reference than the mixed-population class-residualization study.

## What are the principal data structures?

| Data | Concrete representation | Used by |
| --- | --- | --- |
| Expression / abundance matrix | pandas `DataFrame`; genes/proteins × samples or samples × features depending on loader | Cohort-specific preprocessing |
| Feature-space point cloud | NumPy `ndarray`, rows = samples, columns = features or principal components | PCA, ripser, GUDHI |
| Persistence diagrams | `ripser(...)["dgms"]`: arrays of birth/death pairs per homology dimension; GUDHI gives persistence pairs | Maximum finite H1 birth-to-death gap |
| Cohort labels | pandas `Series` or one-dimensional Boolean NumPy array | Tumor/normal comparisons, class residualization |
| Sweep configuration | JSON objects containing tag, imputation rule, feature count, component count, and null-draw counts | `02_run_sweep.py` |
| Monte Carlo results | NumPy arrays, pickled dictionaries, or JSON lists, depending on phase | Null summaries and audit records |
| Audit shards | Compressed NumPy `.npz` with draw indexes and provenance hashes | Clean-lineage worker and aggregator |
| Recorded evidence | CSV, JSON, markdown reports, figures, and `checksums.sha256` | Review and comparison |

These formats are internal research files, **not** a stable package or network API contract. Some historical pickles are executable deserialization inputs; only load copies whose provenance you trust.

## Why is the source organized by cohort rather than one universal pipeline?

The code has modality-specific transforms and exceptions. GSE146889 uses gene-level RNA sequencing data; CPTAC uses log2 protein ratios and per-sample median-centering; TCGA-LUAD methylation clips beta values before a logit transform, uses PCA(35) because of its sample count, and applies a spectral distance for its counted PCA comparison. A single assumed Euclidean PCA(50) path would change the experiment. The cost of cohort-specific modules is duplicated download and preprocessing code, different output filenames, and the need to review each pipeline on its own terms.

## How do validation and provenance work?

`.github/workflows/ci.yml` installs `requirements.txt` with Python 3.11, checks registered SHA-256 hashes, runs the default figure path, checks the independent GSE146889 reconciliation, and runs `pip-audit`. The CPTAC sensitivity workflow separately tests synthetic invariants and runs a reduced-draw pilot on pull requests; 499-draw repetition requires explicit workflow dispatch. Historical procedures are kept under `audit/`, `experiments/`, and `results/` because replacing them would erase the distinction between original measurements and later checks. Refer to [USAGE.md](USAGE.md) for exact commands and [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md) for their costs.
