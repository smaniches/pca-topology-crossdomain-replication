# GSE146889 confound-attribution clean-lineage rerun protocol

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)  
**Status:** frozen before execution  
**Purpose:** independently re-execute the GSE146889 confound-attribution analysis from the public raw GEO matrix without using the transcript-reconstructed driver or its intermediate artifacts.

This is a provenance repair and robustness audit. It does not modify the original cross-domain preregistration or retroactively convert the confound audit into a preregistered study. The Git commit containing this file is the execution-time protocol anchor.

## 1. Independence boundary

The rerun implementation may use the public statistical specification in the committed report and standard published/library algorithms. It must not:

- import `code/confound_attribution_audit/confound_audit_gse146889.py`;
- load `data_gse146889_reused_null/*.pkl`;
- read historical result values until all new statistics have been computed;
- tune preprocessing, seeds, null models, or thresholds to increase agreement with the historical report;
- overwrite or delete historical artifacts.

Historical CSVs are read only after computation to produce a transparent comparison table.

## 2. Data

- Dataset: NCBI GEO `GSE146889`.
- Raw source: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/suppl/GSE146889_GeneCount.tsv.gz`.
- The workflow records the downloaded byte length and SHA-256 digest.
- RPKM columns are those ending in `_rpkm`.
- Rows containing any non-numeric or missing RPKM value are excluded and counted.
- Labels are derived only from `_tumor_` and `_normal_` in the sample names; any ambiguous sample aborts the run.
- Expected cohort cardinality is checked structurally: 176 samples, 91 tumor, 85 normal. This check is not a numerical-outcome check.

## 3. Software and numerical environment

- GitHub-hosted Ubuntu runner.
- Python 3.11.
- Exact package versions from the repository `requirements.txt`.
- The workflow records `python --version`, `pip freeze`, platform metadata, Git commit SHA, script SHA-256, and protocol SHA-256.
- Randomness is generated with NumPy `SeedSequence`/`default_rng` or explicitly recorded `RandomState` seeds as specified below.

## 4. Preprocessing

For the mixed cohort and every independently re-fit subset:

1. `log1p(RPKM)`.
2. Remove zero-variance genes within that analysis subset.
3. Rank by population variance (`ddof=0`) using a stable sort; retain the top 2,000 genes.
4. Record the variance cutoff and the number of genes tied at the cutoff. If a cutoff tie exists, run a documented legacy-order sensitivity calculation; do not choose the result that better matches history.
5. Standardize each retained gene to zero mean and unit variance using `StandardScaler` fit within the analysis subset.
6. Fit `PCA(n_components=min(50, n_samples-1, n_features), random_state=42)` under the pinned scikit-learn version.
7. Compute Vietoris-Rips persistent homology with `ripser`, `maxdim=1`; the primary statistic is the maximum finite H1 death-minus-birth persistence.

## 5. Controls and draw counts

### 5.1 Mixed-set reference

- Observed mixed-cohort PCA(50) maximum H1 persistence.
- Gaussian null: 500 draws, each `N(0,1)` with shape 176 x 2000, standardized and passed through PCA(50) and PH.
- Pipeline-symmetric permutation null: 2,000 draws, independently permuting each selected standardized gene across samples, then re-fitting PCA(50) and PH.
- Base seed: 42.

### 5.2 Within-class decomposition

Re-fit the entire preprocessing pipeline separately for tumor-only and normal-only samples.

For each class:

- Gaussian null: 300 draws, matched sample count x 2000, standardized, PCA(50), PH.
- Permutation null: 300 draws on the class-specific selected standardized genes, PCA re-fit each draw, PH.
- Seeds: tumor Gaussian 14688901; normal Gaussian 14688902; tumor permutation 14688911; normal permutation 14688912.

### 5.3 Within-stratum control

- Compute the unit class-mean-shift direction in the mixed standardized HVG space.
- Bin the continuous projection into four equal-count quartiles using deterministic rank order.
- Re-fit the full preprocessing pipeline separately within each quartile.
- For each quartile: 300 Gaussian and 300 permutation draws, matched to the quartile sample count.
- Seeds: Gaussian `14689000 + quartile_index`; permutation `14689100 + quartile_index`.

### 5.4 Residualization

- In mixed standardized HVG space, shift each sample by subtracting its class-specific mean and adding the pooled grand mean, exactly equalizing class means gene by gene.
- Re-fit PCA(50) and PH on the residualized matrix.
- Compare intact and residualized statistics to the newly generated mixed-set Gaussian and permutation nulls from Section 5.1.
- Compute five-fold stratified CV AUC under two explicitly reported logistic-regression regularizations: `C=0.01` (historical confound-audit implementation) and `C=1.0` (sensitivity).
- Run a two-sided label-permutation test on residualized-space AUC with 200 permutations, seed 14689300.
- Request H1 cocycles from `ripser`; for the most persistent H1 class, report the support vertices of the corresponding cocycle and their tumor/normal composition. This is labelled a cocycle-support diagnostic, not a canonical homology-cycle identity.

### 5.5 Block bootstrap

For mixed, tumor-only, and normal-only PCA point clouds:

- 2,000 row-bootstrap resamples with replacement.
- Recompute maximum H1 persistence on each resample.
- Report percentile 95% intervals for `observed - matched Gaussian-null mean` and its z-scaled version.
- Seeds: mixed 14689200; tumor 14689201; normal 14689202.
- The Gaussian reference distribution is held fixed, matching the historical audit; this approximation remains explicitly disclosed.

## 6. Outputs

The run must write all computed quantities before historical comparison:

- `run_provenance.json`;
- `primary_results.json`;
- `within_class.csv`;
- `within_stratum.csv`;
- `residualization.csv`;
- `bootstrap.csv`;
- `confound_spectrum.csv`;
- `null_and_bootstrap_draws.npz`;
- `gse146889_clean_lineage_audit.png`, generated only from the new machine-readable outputs.

After those objects exist in memory, the script may read the historical committed CSVs and write:

- `historical_comparison.csv`;
- `REPORT.md`, preserving all disagreements.

## 7. Interpretation rules

The rerun is not required to match every historical Monte Carlo value exactly. Assessment distinguishes:

- deterministic real-data agreement;
- Monte Carlo agreement within sampling variability;
- qualitative conclusion agreement;
- material disagreement that changes a scientific claim.

The historical qualitative verdict is sustained only if the clean rerun supports both of these statements:

1. the aggregate topological excess is not eliminated by linear class-mean residualization and persists in most confound-projection strata; and
2. single-class robustness is weaker than the mixed result, including bootstrap uncertainty that may cross zero.

If either statement fails, the report must narrow or reverse the historical verdict. No result is forced to agree.

## 8. Execution and publication boundary

The workflow may commit generated outputs only to `audit/gse146889-clean-lineage`. It must not merge, tag, release, alter the manuscript, or modify `main`. Integration requires a separate exact-head review and explicit authorization.