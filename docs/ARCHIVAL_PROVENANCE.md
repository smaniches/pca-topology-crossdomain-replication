# GitHub–Zenodo archive provenance

**Updated verification:** 2026-10-11 01:30 UTC. The v0.3.0 GitHub tag, Zenodo version and shared concept DOI, downloaded ZIP, all registered source-file SHA-256 hashes, the matched 499-draw numerical record and original PDFs have passed independent read-only GitHub Actions verification. This page distinguishes original research records from later corrective evidence and provides the authoritative DOI for each version.

## Which releases and DOIs correspond?

| Object | Identifier | Verified metadata |
| --- | --- | --- |
| Original GitHub release | [`v0.1.0-biorxiv`](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.1.0-biorxiv) | Published 2026-07-10, historical source commit `4fa4ddee67e529e170cd031f3852a6875b40d9e5` |
| Original Zenodo version | [**10.5281/zenodo.21287944**](https://doi.org/10.5281/zenodo.21287944) | Zenodo record `21287944`, version `v0.1.0-biorxiv`, linked to this GitHub repository |
| Corrected GitHub release | [**v0.2.0**](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0) | Published **2026-10-10T00:54:55Z** (October 9, Eastern time), source commit `c2f227884c360b6d7ce9db4056c4a0ff14a87f1f`, GitHub release ID `408470234` |
| v0.2.0 Zenodo version | [**10.5281/zenodo.23272317**](https://doi.org/10.5281/zenodo.23272317) | Zenodo record `23272317`, version `v0.2.0`; immutable prior corrected snapshot |
| Finalized matched-study GitHub release | [**v0.3.0**](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.3.0) | Published **2026-10-11T01:28:23Z** (October 10 in US Eastern time), source commit `bb453a01874a97d996bdd19a7eb571c3c2bc82ae`, GitHub release ID `409320182` |
| v0.3.0 Zenodo version | [**10.5281/zenodo.23290337**](https://doi.org/10.5281/zenodo.23290337) | **Verified** new Zenodo record `23290337`, version `v0.3.0`, matched source association and downloaded ZIP SHA-256 registry |
| **All-versions concept DOI** | [**10.5281/zenodo.21287943**](https://doi.org/10.5281/zenodo.21287943) | All three published Zenodo versions report **`conceptrecid = 21287943`**. The v0.3.0 study is a new version within the same family, **not** an unrelated deposit |

The [v0.2.0 release notes](../release/v0.2.0/RELEASE_NOTES.md) and [v0.3.0 matched-study release notes](../release/v0.3.0/RELEASE_NOTES.md) record the distinct scientific scopes. Version-specific Zenodo DOIs were discovered and verified **after** their respective GitHub releases. Source tags are immutable; an archive's tagged `CITATION.cff` may predate DOI minting. The current `main` branch's [`CITATION.cff`](../CITATION.cff) now cites the subsequently verified **v0.3.0 version DOI**. The original July and corrected v0.2 version identifiers remain unchanged.

## Did automatic archiving really happen?

**Yes.** The repository owner's GitHub–Zenodo integration ingested the new [GitHub Release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0). The [public Zenodo verification workflow](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38011308564) obtained both records from Zenodo's live API and repository DOI badge. Its archived machine-readable report identifies the exact record IDs, concept family, versions, public file checksums, and filenames.

Zenodo returned the following earlier-version ZIP metadata (v0.3.0 ZIP evidence is documented below):

| Version | Archive filename | Size (bytes) | Zenodo MD5 |
| --- | --- | ---: | --- |
| July `v0.1.0-biorxiv` | `smaniches/pca-topology-crossdomain-replication-v0.1.0-biorxiv.zip` | 11,260,952 | `1200f58533ed56b46b310f832e0dfc95` |
| October `v0.2.0` | `smaniches/pca-topology-crossdomain-replication-v0.2.0.zip` | 11,654,197 | `d3dcb72a6fd9481cb1e1a137b90e160a` |

The archived ZIP byte-level comparison is implemented separately by [`audit/verify_zenodo_archive_bytes.py`](../audit/verify_zenodo_archive_bytes.py). **Independent archive-byte verification passed.** The read-only [verification workflow](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38011516246) downloaded both Zenodo ZIPs, checked each file against its Zenodo-declared MD5, recomputed the ZIP SHA-256, and checked the original sealed PDFs. For the v0.2.0 archive it also verified **all 188 registered file SHA-256 entries** against the ZIP's actual source files (189 files total including the registry). The result is saved in [GitHub Actions artifact 11653925732](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38011516246/artifacts/11653925732). Verified archive SHA-256 digests: July `bafb3bf5eb0c3cc96d966a223f2f4a129514468808fa1700a31e410e2d2d39e4`; October `6f58a6bbed54a247574134adb901d618398dd7ca656afccacbc0cb8ccbf7c9c7`.

Zenodo integrates **GitHub Releases**, not every commit to `main`. Subsequent documentation updates are not automatically present in the frozen v0.2.0 Zenodo ZIP. The 2026-10-09 release snapshot, not the later `main` head, is the authoritative byte set for this deposit.

## How should each body of evidence be cited?

For the **original preregistered work**, cite DOI [10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944) and its historical July GitHub release when source provenance matters. The original [preregistration](../prereg/PREREGISTRATION.md), submission package, submitted PDFs, and historical result files are preserved.

For the **prior v0.2.0 corrected research-software snapshot and its October 9 experiments**, cite DOI **[10.5281/zenodo.23272317](https://doi.org/10.5281/zenodo.23272317)** and its [v0.2.0 GitHub source tag](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0). This release includes the [methodological correction](../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md), the [CPTAC mixed-cohort sensitivity](../results/cptac_resid_null_sensitivity_20261009/REPORT.md), the [tumor-only covariance-preserving sensitivity](../results/cptac_tumor_covariance_null_20261009/REPORT.md), and [all 499 tumor-only null draws](../results/cptac_tumor_covariance_null_20261009/499_draws.json). For stable citation to the **whole series**, use the verified concept DOI [10.5281/zenodo.21287943](https://doi.org/10.5281/zenodo.21287943), while naming the specific version when reproducing numbers.

Scientific qualification remains essential: the historical classifier area-under-the-curve permutation `p=0.005` claims are withdrawn, the methylation result is metric-gated, and the new tumor-only covariance-preserving experiment did **not** reject its null (`p=0.092`). The earlier mixed 194-sample `p=0.002` control is not a matched null-model comparison with that 110-sample tumor-only experiment.

## Closed v0.3.0 matched comparison — verified archive

The fully executed and independently replayed [110-tumor matched study](../results/cptac_tumor_matched_null_20261010/REPORT.md) and [all 499 new feature-permutation draws](../results/cptac_tumor_matched_null_20261010/499_draws.json) are now published in the **[GitHub v0.3.0 release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.3.0)** and its independently verified **[Zenodo version DOI 10.5281/zenodo.23290337](https://doi.org/10.5281/zenodo.23290337)**. See the [v0.3.0 release notes](../release/v0.3.0/RELEASE_NOTES.md). The new version is linked to the existing [concept DOI 10.5281/zenodo.21287943](https://doi.org/10.5281/zenodo.21287943) rather than creating an unrelated family.

**Exact scientific outcome:** for the same 110 CPTAC tumor samples, 2,000-protein preprocessed panel, PCA(50), Euclidean persistent H1 and observed value 5.3542976379, 0/499 independent protein-column permutations met/exceeded the statistic (**p = 0.002**), versus 45/499 full-covariance-preserving Haar controls (**p = 0.092**). This is a **post-preregistration, outcome-informed** matched reference-model sensitivity analysis; it does **not** establish biological topology or exclusive covariance causation.

**Independent verification:** [GitHub Actions run 38102031352](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38102031352) checked new Zenodo record `23290337` (same concept family), downloaded the v0.3.0 ZIP (**11,695,814 bytes**), verified its Zenodo-declared MD5, independently recomputed its ZIP SHA-256, validated **all 200 tracked-file SHA-256 entries** inside the immutable source archive, and checked the original sealed PDF digests and complete archived 499-draw evidence. The [machine-readable verification artifact](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38102031352) preserves full details. The earlier v0.2.0 ZIP and DOI were also checked and left unchanged.

**Immutable snapshot boundary:** later `main` documentation/CFF changes containing the minted DOI are **not automatically inside the v0.3.0 GitHub tag or Zenodo ZIP**, because those source archives were frozen at publication. The post-release DOI citation is a later metadata correction, not a reason to overwrite the validated archive.

The separately documented [independent-cohort validation phase](INDEPENDENT_COHORT_VALIDATION_PHASE.md) remains **NOT STARTED**. Publishing the matched null experiment does not constitute, initiate, or complete independent biological replication.

## Historical submission and legal boundaries

The original bioRxiv submission `BIORXIV/2026/737609` was **not posted**, following an organizational-affiliation screening rule. This was not a scientific peer-review rejection. Statements within `submission/biorxiv/` saying that a DOI had not yet been issued are historical statements from package creation and must not be rewritten.

The [code license](../LICENSE) is MIT for code; the original manuscript text and figures are expressly excluded by its license note. Zenodo record metadata cannot override the source material's actual rights without independent authority.

## Preservation rules

- Never delete, overwrite, retag, or retroactively change either GitHub release or an already published Zenodo ZIP.
- Do not use a version-specific DOI to identify a different release.
- After each **new GitHub Release**, independently verify Zenodo's version DOI, concept DOI, source association, and archived file integrity. Do not infer synchronization from a successful GitHub CI run alone.
- Keep correction notices available alongside any references to the original manuscript.
- The read-only [Zenodo record probe](../audit/probe_zenodo_v0.2.0.py) and [archive-byte verifier](../audit/verify_zenodo_archive_bytes.py) can be re-executed without creating a deposit.

**Archive state: VERIFIED.** Both Zenodo records, their shared concept DOI, the original release's sealed PDFs, the actual archived ZIP bytes, and all registered v0.2.0 source-file checksums passed independent public-network verification. This does not imply that future `main` commits automatically appear in the frozen v0.2.0 archive.
