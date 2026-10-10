# GitHub–Zenodo archive provenance

**Verified:** 2026-10-10 00:58 UTC by a public read-only Zenodo API and GitHub Actions probe. This page distinguishes original research records from later corrective evidence and provides the authoritative DOI for each version.

## Which releases and DOIs correspond?

| Object | Identifier | Verified metadata |
| --- | --- | --- |
| Original GitHub release | [`v0.1.0-biorxiv`](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.1.0-biorxiv) | Published 2026-07-10, historical source commit `4fa4ddee67e529e170cd031f3852a6875b40d9e5` |
| Original Zenodo version | [**10.5281/zenodo.21287944**](https://doi.org/10.5281/zenodo.21287944) | Zenodo record `21287944`, version `v0.1.0-biorxiv`, linked to this GitHub repository |
| Corrected GitHub release | [**v0.2.0**](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0) | Published **2026-10-10T00:54:55Z** (October 9, Eastern time), source commit `c2f227884c360b6d7ce9db4056c4a0ff14a87f1f`, GitHub release ID `408470234` |
| New Zenodo version | [**10.5281/zenodo.23272317**](https://doi.org/10.5281/zenodo.23272317) | Zenodo record `23272317`, version `v0.2.0`, associated with this repository by its public record metadata |
| **All-versions concept DOI** | [**10.5281/zenodo.21287943**](https://doi.org/10.5281/zenodo.21287943) | Both Zenodo versions report **`conceptrecid = 21287943`**. This is the same version family, **not** two unrelated deposits |

The new Zenodo DOI was discovered and verified **after** the GitHub release was published. The original July DOI has not been repurposed. The v0.2.0 source tag is immutable and contains the citation file *as it existed before Zenodo minted the new DOI*. The current `main` branch's [`CITATION.cff`](../CITATION.cff) now correctly cites the subsequently verified **v0.2.0 version DOI**, so it differs from the earlier tagged source in this expected post-publication metadata field.

## Did automatic archiving really happen?

**Yes.** The repository owner's GitHub–Zenodo integration ingested the new [GitHub Release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0). The [public Zenodo verification workflow](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38011308564) obtained both records from Zenodo's live API and repository DOI badge. Its archived machine-readable report identifies the exact record IDs, concept family, versions, public file checksums, and filenames.

Zenodo returned these ZIP file metadata:

| Version | Archive filename | Size (bytes) | Zenodo MD5 |
| --- | --- | ---: | --- |
| July `v0.1.0-biorxiv` | `smaniches/pca-topology-crossdomain-replication-v0.1.0-biorxiv.zip` | 11,260,952 | `1200f58533ed56b46b310f832e0dfc95` |
| October `v0.2.0` | `smaniches/pca-topology-crossdomain-replication-v0.2.0.zip` | 11,654,197 | `d3dcb72a6fd9481cb1e1a137b90e160a` |

The archived ZIP byte-level comparison is implemented separately by [`audit/verify_zenodo_archive_bytes.py`](../audit/verify_zenodo_archive_bytes.py). **A returned Zenodo checksum is verified metadata, not proof of downloaded-byte integrity until the corresponding download and independent hash verification pass.** The verification [workflow](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/workflows/verify-zenodo-v0.2.0.yml) records that distinction and saves the results as artifacts.

Zenodo integrates **GitHub Releases**, not every commit to `main`. Subsequent documentation updates are not automatically present in the frozen v0.2.0 Zenodo ZIP. The 2026-10-09 release snapshot, not the later `main` head, is the authoritative byte set for this deposit.

## How should each body of evidence be cited?

For the **original preregistered work**, cite DOI [10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944) and its historical July GitHub release when source provenance matters. The original [preregistration](../prereg/PREREGISTRATION.md), submission package, submitted PDFs, and historical result files are preserved.

For the **current corrected research-software snapshot and October experiments**, cite DOI **[10.5281/zenodo.23272317](https://doi.org/10.5281/zenodo.23272317)** and its [v0.2.0 GitHub source tag](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0). This release includes the [methodological correction](../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md), the [CPTAC mixed-cohort sensitivity](../results/cptac_resid_null_sensitivity_20261009/REPORT.md), the [tumor-only covariance-preserving sensitivity](../results/cptac_tumor_covariance_null_20261009/REPORT.md), and [all 499 tumor-only null draws](../results/cptac_tumor_covariance_null_20261009/499_draws.json). For stable citation to the **whole series**, use the verified concept DOI [10.5281/zenodo.21287943](https://doi.org/10.5281/zenodo.21287943), while naming the specific version when reproducing numbers.

Scientific qualification remains essential: the historical classifier area-under-the-curve permutation `p=0.005` claims are withdrawn, the methylation result is metric-gated, and the new tumor-only covariance-preserving experiment did **not** reject its null (`p=0.092`). The earlier mixed 194-sample `p=0.002` control is not a matched null-model comparison with that 110-sample tumor-only experiment.

## Historical submission and legal boundaries

The original bioRxiv submission `BIORXIV/2026/737609` was **not posted**, following an organizational-affiliation screening rule. This was not a scientific peer-review rejection. Statements within `submission/biorxiv/` saying that a DOI had not yet been issued are historical statements from package creation and must not be rewritten.

The [code license](../LICENSE) is MIT for code; the original manuscript text and figures are expressly excluded by its license note. Zenodo record metadata cannot override the source material's actual rights without independent authority.

## Preservation rules

- Never delete, overwrite, retag, or retroactively change either GitHub release or an already published Zenodo ZIP.
- Do not use a version-specific DOI to identify a different release.
- After each **new GitHub Release**, independently verify Zenodo's version DOI, concept DOI, source association, and archived file integrity. Do not infer synchronization from a successful GitHub CI run alone.
- Keep correction notices available alongside any references to the original manuscript.
- The read-only [Zenodo record probe](../audit/probe_zenodo_v0.2.0.py) and [archive-byte verifier](../audit/verify_zenodo_archive_bytes.py) can be re-executed without creating a deposit.

**Archive state:** Both v0.1.0 and v0.2.0 Zenodo records are **verified**, and their shared DOI family is **verified** by direct public metadata. File byte-level integrity is a separate check with independently recorded evidence.
