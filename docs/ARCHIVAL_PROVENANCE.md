# Archive provenance and citation status

**Reviewed:** 2026-10-09. **Scope:** the public GitHub repository `smaniches/pca-topology-crossdomain-replication`. This record tracks what can and cannot be established about the published archive. It is not a new Zenodo deposit.

## Which versions exist in the available evidence?

| Object | Dated evidence | Safe interpretation |
| --- | --- | --- |
| Locked original protocol | [Preregistration](../prereg/PREREGISTRATION.md), July 2026 | Original hypotheses and decisions remain unchanged |
| Original GitHub release | [`v0.1.0-biorxiv`](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.1.0-biorxiv), **2026-07-10**, release commit `4fa4ddee67e529e170cd031f3852a6875b40d9e5` | Preserved historical GitHub snapshot, prior to October corrections |
| Original bioRxiv submission | `BIORXIV/2026/737609`; [submission seal](../submission/biorxiv/SUBMISSION_SEAL.md) | The manuscript was **not published on bioRxiv**. The 2026-07-14 screening response cited lack of an accepted organizational affiliation, **not an assessment of scientific merit** |
| DOI recorded in GitHub | [`CITATION.cff`](../CITATION.cff) and [2026-08-20 citation PR #4](https://github.com/smaniches/pca-topology-crossdomain-replication/pull/4) | The repository records **[10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944)** as a prior citable reference. This alone does **not** verify the current Zenodo record, its files, or its DOI version type |
| October corrections and sensitivity experiments | [Methodological corrections](../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md); [CPTAC experiment](../results/cptac_tumor_covariance_null_20261009/REPORT.md); [numeric evidence](../results/cptac_tumor_covariance_null_20261009/499_draws.json) | New GitHub research added after the July release. Its inclusion in any Zenodo record has **not** been established |

The original `submission/biorxiv/` files contain statements such as “no Zenodo DOI exists yet.” Those statements were made **at package-creation time**, before the later release and citation change. They are preserved as historical records; editing them to describe today's state would damage their provenance.

## Does Zenodo automatically mirror current GitHub files?

**No.** Zenodo's optional GitHub integration, when enabled, archives **new GitHub releases**, not each commit on `main`. A changed README, a merged pull request, and a newly committed scientific result do **not** prove that any existing Zenodo deposit has changed. See the official [GitHub integration guidance](https://help.zenodo.org/docs/github/enable-repository/).

A Zenodo **version DOI** identifies a specific archived version. A **concept DOI** refers to a version family and can resolve to its latest version. Metadata-only edits generally do not require a new version, but replacing or adding archived files does. See the official [DOI versioning guidance](https://help.zenodo.org/docs/deposit/manage-versions/).

**Unverified for this project:** the live DOI landing page and files, its publication date, the DOI's concept-versus-version status, all other versions or independently made mirrors, Zenodo account ownership, file hashes, and whether the GitHub integration is enabled. No direct Zenodo/API record could be retrieved during this review. This repository exposes **one DOI in its tracked citation metadata**, but that cannot establish how many separate Zenodo records exist.

## How should the project be cited now?

For the **original archived work**, start with DOI [10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944), as recorded by the project, and **verify its actual Zenodo landing page, title, files, and version** before relying on it. The historical GitHub source can independently be cited by release tag `v0.1.0-biorxiv` and its commit SHA.

For the **October corrections**, cite [the correction notice](../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md) **alongside** the original work. The correction withdraws unsupported residualized-classifier permutation significance claims and clarifies that the methylation PCA(35) zero result was obtained under a changed spectral distance.

For the **October CPTAC covariance-preserving result**, cite the [frozen protocol](../audit/cptac_covariance_null/PROTOCOL.md), [full report](../results/cptac_tumor_covariance_null_20261009/REPORT.md), [499 raw values](../results/cptac_tumor_covariance_null_20261009/499_draws.json), and the Git commit containing them. **Do not cite the earlier DOI as though it already archived this October experiment.** A new DOI or linked record must be verified after publication.

## What needs to happen before updating Zenodo?

1. In the **owner's Zenodo account**, open [the recorded DOI](https://doi.org/10.5281/zenodo.21287944). Identify the concept DOI and every version-specific DOI. Search for any **other manually created mirrors** before creating a new deposit.
2. Inspect and record each version's filename list, checksums, publication date, title, release tag or linked GitHub commit, and licensing. Compare archived files against the original GitHub release; do not assume equality.
3. If a correction link in existing metadata is appropriate, add a dated description or related identifier using the account's record-editing functionality. **Do not replace the original sealed PDF or old files.**
4. To archive the newer corrected research files, first confirm whether a **new version of the existing Zenodo version family** or a separate intentionally related record is appropriate. Preserve the original version DOI. Do not trigger a new GitHub release until the integration and target-version relationship are confirmed; otherwise a second, unrelated Zenodo family may be created.
5. After any authorized publication, confirm the resulting version DOI, concept DOI, visible files, hashes and cross-links. Only then update [`CITATION.cff`](../CITATION.cff) and the project citation instructions with the exact verified identifiers.

**Change-control status:** This documentation update does **not** edit Zenodo metadata or records, create a Zenodo deposit, publish a GitHub release or tag, change the preregistration, modify the submitted PDF, or recompute historical scientific results. The external Zenodo content/versions remain **UNKNOWN** until directly inspected.
