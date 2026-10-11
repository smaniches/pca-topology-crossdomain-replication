# Principal component analysis (PCA) and persistent homology (PH): cross-cohort replication

## What this is and what problem it solves

This is a research codebase for checking whether a change in the longest-lived one-dimensional homology (H1) feature after PCA reflects the data or also appears in generated controls. It contains a pilot analysis, cohort-specific replication scripts, a parameter sweep, confound checks, and independently recorded follow-up audits. The command-line entry point is `reproduce.py`; it is not a hosted service or a reusable inference library.

**Publication and citation:** The original [July 2026 GitHub release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.1.0-biorxiv) is preserved under Zenodo version DOI [10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944). The corrected October research was released on GitHub as [**v0.2.0**](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0), source commit `c2f227884c360b6d7ce9db4056c4a0ff14a87f1f`, and **automatically archived by Zenodo as [10.5281/zenodo.23272317](https://doi.org/10.5281/zenodo.23272317)**. Both belong to the same version family, [concept DOI 10.5281/zenodo.21287943](https://doi.org/10.5281/zenodo.21287943). The original bioRxiv manuscript was submitted but **not posted** following affiliation screening. See [v0.2.0 release notes](release/v0.2.0/RELEASE_NOTES.md), the [new v0.3.0 matched-experiment release notes](release/v0.3.0/RELEASE_NOTES.md), and [archive provenance and verification](docs/ARCHIVAL_PROVENANCE.md). The completed matched study is **published separately as [GitHub v0.3.0](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.3.0)** and **verified at Zenodo version DOI [10.5281/zenodo.23290337](https://doi.org/10.5281/zenodo.23290337)**, which belongs to the same [all-versions concept DOI family](https://doi.org/10.5281/zenodo.21287943). Its archived source commit is `bb453a01874a97d996bdd19a7eb571c3c2bc82ae`.

**Latest matched scientific result (2026-10-10):** On the **identical 110 tumor samples, 2,000 features, PCA(50), and H1 statistic**, the observed value **5.354** exceeded all 499 independently shuffled-protein controls (**p = 0.002**), but not the 499 controls preserving complete cross-protein covariance at the 0.05 threshold (**45 exceedances; p = 0.092**). This supports dependence-sensitive geometry, **not** an established biological topology beyond covariance. The two nulls preserve different properties, so their p-values alone do not prove covariance causation. [Read the fully executed comparison](results/cptac_tumor_matched_null_20261010/REPORT.md) and [all 499 new values](results/cptac_tumor_matched_null_20261010/499_draws.json). This study is **included in the independently verified v0.3.0 Zenodo ZIP** and is not part of the earlier immutable v0.2.0 snapshot. **Independent-cohort replication is future research**, not unfinished work on this comparison; see the [separate phase plan](docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md).

![Actual histogram of the covariance-preserving reference, one arm of the two-null comparison.](results/cptac_tumor_covariance_null_20261009/null_distribution.svg)

*This histogram shows the original covariance-preserving arm only; the newly matched independent-feature reference and its full distribution are documented in the comparison report.*

## Why existing approaches fall short

Comparing the observed H1 persistence before and after PCA alone cannot distinguish a data-specific change from one induced by the projection. The scripts therefore compare the observed statistic with Gaussian and column-permutation controls processed through corresponding analysis steps. Those controls have their own assumptions; agreement or disagreement does not identify biological mechanisms.

## What is genuinely new and what is inherited

The project-specific contribution is a locked, cross-cohort comparison protocol and its executable reproduction and audit records, including a separate class-conditional sensitivity analysis of Clinical Proteomic Tumor Analysis Consortium (CPTAC) data. The implementation inherits PCA and classification from scikit-learn, numerical processing from NumPy and pandas, persistent-homology computation from ripser and GUDHI, and public cohort data from the source archives identified in the scripts. It does not introduce a new PCA algorithm, homology algorithm, or proof that observed loops are biological. The original protocol is at [prereg/PREREGISTRATION.md](prereg/PREREGISTRATION.md); later experiments are labeled separately.

## Quickstart

On Linux, install Git, Python 3.11 with `venv`, and `sha256sum`. Network access is needed to clone and install packages, but the default run uses committed pilot checkpoints and needs no data-service credentials. From a Bash shell, execute these commands in order:

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

A successful run prints `REPRODUCE: SUCCESS` and leaves two nonempty image files under `reproduce_output/`. This checks the default figure-generation path, **not** the raw-data analyses or the paper's claims.

## Where it breaks

The default path unpickles committed checkpoints and must not be used with untrusted replacements. Full-data paths depend on external Gene Expression Omnibus, Genomic Data Commons, or Proteomic Data Commons services and can fail when data or metadata change. Reducing null draws changes reported tail probabilities. The methylation branch changes its distance metric, and historical label-conditioned residualization is not a valid held-out predictive evaluation; the later CPTAC sensitivity analysis fixes only one distinct null-comparison question. See the documented failure modes before using results as evidence. **Before citing the existing manuscript or conclusions, read the [2026-10-09 methodological correction](paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md).**

## Documentation

- [Architecture and data flow](docs/ARCHITECTURE.md)
- [Design choices and costs](docs/DESIGN_DECISIONS.md)
- [Limits and known failures](docs/LIMITATIONS.md)
- [Commands, data sources, outputs, and troubleshooting](docs/USAGE.md)
- [Publication status, DOI versions, and Zenodo provenance](docs/ARCHIVAL_PROVENANCE.md)
