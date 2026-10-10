# v0.2.0 — Methodological corrections and reproducible follow-up research

**Released:** 2026-10-09  
**Scope:** a versioned research-software and evidence update to the earlier `v0.1.0-biorxiv` GitHub release. The July submission package and its submitted PDF remain historical, byte-preserved artifacts.

## Why this version exists

The original cross-domain persistent-homology project needed more precise claims about its statistical controls. This version **does not rewrite** the original preregistration or historical results. Instead, it makes the identified methodological limitations explicit, adds independently executed sensitivity analyses, and supplies reproducible records and tested onboarding documentation.

## Scientific corrections that change interpretation

- **Residualized classifier AUC significance is withdrawn.** The historical class-mean residualization uses the true labels of all observations before cross-validation. Its label-permutation test does not recompute the label-conditioned features on each permutation. Previously reported `p=0.005` significance values from that procedure are not valid evidence of independent confound removal; the historical numerical AUC records remain unchanged.
- **Historical post-residualization nulls were not transform matched.** Class-mean centering was applied to real observations but not identically to every Gaussian/feature-permutation surrogate.
- **TCGA-LUAD methylation depends on the spectral distance.** On 36 samples, PCA(35) preserves centered Euclidean pairwise distances. Its zero H1 observation occurs under the different preregistered spectral distance, so it is not evidence that Euclidean PCA itself destroyed a loop.

See [the full methodological correction](../../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md).

## Executed follow-up experiments

1. **CPTAC mixed-cohort conditional feature-permutation sensitivity.** All 499 null values are archived. The observed class-residualized PCA(50) maximum H1 persistence was 4.40994, with zero of 499 shuffled-feature nulls reaching it; one-sided plus-one `p=0.002`. This reference destroys cross-protein covariance and does not validate biological attribution. [Report](../../results/cptac_resid_null_sensitivity_20261009/REPORT.md) · [499 draws](../../results/cptac_resid_null_sensitivity_20261009/499_draws.json).
2. **CPTAC tumor-only covariance-preserving sensitivity.** On 110 tumor samples, all 499 Haar-rotation null draws preserve the full sample-centered feature covariance matrix. Observed maximum H1 was 5.35430; 45 of 499 null values reached or exceeded it; plus-one `p=0.092`. **This does not meet the prespecified 0.05 sensitivity threshold.** [Frozen protocol](../../audit/cptac_covariance_null/PROTOCOL.md) · [Report](../../results/cptac_tumor_covariance_null_20261009/REPORT.md) · [499 draws](../../results/cptac_tumor_covariance_null_20261009/499_draws.json) · [data-derived null histogram](../../results/cptac_tumor_covariance_null_20261009/null_distribution.svg).

**These two p-values must not be contrasted as a controlled experiment.** They use different populations (194 mixed tumor/normal after residualization versus 110 tumor-only) and different null constructions. No biological-topology mechanism is established by this release.

## Documentation, provenance, and checks

- A shorter [README](../../README.md) with a clean Python 3.11 virtual-environment quickstart; full [architecture](../../docs/ARCHITECTURE.md), [design rationale](../../docs/DESIGN_DECISIONS.md), [limitations](../../docs/LIMITATIONS.md), and [usage guide](../../docs/USAGE.md).
- [Plain-language experiment explanation](../../docs/EXPERIMENT_EXPLAINED.md), plus an SVG chart derived deterministically from the entire 499-draw covariance-preserving reference distribution.
- [GSE146889 clean-lineage reconciliation](../../audit/gse146889_clean_lineage/RECONCILIATION.md), separately executable CPTAC full-count replay, archived raw null values, and immutable original submission-PDF hashes.
- [Archive-provenance guide](../../docs/ARCHIVAL_PROVENANCE.md) and [citation metadata](../../CITATION.cff) distinguish the original July release from October updates.
- GitHub Actions verifies SHA-256 registered files, the original sealed PDFs, quickstart outputs, the data-derived histogram, GSE146889 reconciliation, and pinned dependency security. Research-specific workflows run unit/invariance checks and bounded null pilots.

## Source and archival identity

- **New GitHub release tag:** `v0.2.0`. It should point to the exact successful CI-tested source commit used for this release; never substitute a floating branch or an older July tag.
- **Original release:** [`v0.1.0-biorxiv`](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.1.0-biorxiv) (2026-07-10). Its original DOI is recorded as [10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944) in the provenance guide. This **is not** the DOI of v0.2.0 unless Zenodo independently says so.
- **New Zenodo archive:** the repository owner has enabled GitHub–Zenodo integration, so GitHub Release publication is expected to trigger Zenodo ingestion. The resulting DOI, record, version family, uploaded files, and successful ingestion **must be verified at Zenodo after publication**. Until verified, cite this version using its GitHub release URL/tag and source commit; do not reuse the July DOI as if it identifies the October snapshot.

## What this release deliberately does not do

It does not make the original bioRxiv submission publicly posted, modify the original preregistration, replace the July archived files, certify any confound independence claim from an invalid AUC test, establish a biological mechanism, or demonstrate whether the observed difference between two nonidentical CPTAC null experiments is caused by covariance alone.

**Historical submission status:** the original bioRxiv submission was not posted because of a screening rule regarding organizational affiliation, not a substantive finding about this analysis.
