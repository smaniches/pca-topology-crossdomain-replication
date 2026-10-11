"""Check immutable submission evidence and live-documentation citation routing.

This is an *offline* repository regression check. It intentionally makes no
claim that the Zenodo service, DOI landing page, or any external mirror is
currently accessible. That requires a separate direct archival inspection.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL_DOI = "10.5281/zenodo.21287944"
PREVIOUS_DOI = "10.5281/zenodo.23272317"
CURRENT_DOI = "10.5281/zenodo.23290337"
CONCEPT_DOI = "10.5281/zenodo.21287943"
ORIGINAL_PDF_SHA256 = "cfc2410130c5872979a62e08cecfb964aa5412d42192c0c632c22d48e1983aac"
ORIGINAL_RELEASE = "v0.1.0-biorxiv"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"ARCHIVE PROVENANCE CHECK FAILED: {message}")


def read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def main() -> None:
    citation = read("CITATION.cff")
    readme = read("README.md")
    ledger = read("docs/ARCHIVAL_PROVENANCE.md")
    correction = read("paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md")
    experiment = read("results/cptac_tumor_covariance_null_20261009/REPORT.md")

    # v0.3.0 represents a distinct immutable source snapshot. Its new Zenodo
    # version DOI is not known before GitHub Release publication.
    require(
        'version: 0.3.0' in citation.splitlines(),
        "current CITATION.cff does not declare version 0.3.0",
    )
    require(
        'date-released: 2026-10-10' in citation.splitlines(),
        "current citation date is not 2026-10-10",
    )
    current_doi_lines = [line for line in citation.splitlines()
                         if line.startswith("doi:")]
    for historical in (ORIGINAL_DOI, PREVIOUS_DOI, CONCEPT_DOI):
        require(not any(historical in line for line in current_doi_lines),
                "a historical Zenodo DOI is incorrectly assigned to v0.3.0")
        require(historical in readme and historical in ledger,
                "a historical DOI is missing from the public provenance notes")
    # A new v0.3.0 version DOI is now independently verified. The tagged
    # source predates that DOI, but current main citation metadata must not.
    require(current_doi_lines == ['doi: "10.5281/zenodo.23290337"'],
            "current CFF does not identify the verified v0.3.0 DOI")
    require(CURRENT_DOI in ledger and CURRENT_DOI in readme,
            "v0.3.0 minted DOI is absent from documentation")
    require("version: 0.3.0" in citation.splitlines(),
            "CFF release version is not 0.3.0")
    require("date-released: 2026-10-10" in citation.splitlines(),
            "CFF release date is incorrect")
    require("release/v0.3.0/RELEASE_NOTES.md" in readme,
            "new version release notes missing from README")
    require("release/v0.3.0/RELEASE_NOTES.md" in ledger,
            "new version release notes missing from archive ledger")
    require("docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md" in readme,
            "separate, unexecuted independent-cohort phase is undocumented")
    require(ORIGINAL_DOI in readme, "original DOI omitted from README")
    require("release/v0.2.0/RELEASE_NOTES.md" in readme,
            "README omits v0.2.0 release notes")
    require("release/v0.2.0/RELEASE_NOTES.md" in ledger,
            "archive provenance omits v0.2.0 release notes")
    require("docs/ARCHIVAL_PROVENANCE.md" in readme, "provenance ledger omitted from README")
    require(ORIGINAL_DOI in ledger, "ledger does not identify existing DOI")
    require(ORIGINAL_RELEASE in ledger, "ledger does not pin historical GitHub release")
    require("paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md" in ledger,
            "ledger does not point to methodological correction")
    require("ARCHIVAL_PROVENANCE.md" in correction,
            "correction is missing its version-provenance pointer")
    require("ARCHIVAL_PROVENANCE.md" in experiment,
            "October experiment is missing its version-provenance pointer")

    # Never silently overwrite the original submitted manuscript while
    # refreshing newer README files, supplementary reports, or hashes.
    for relative in ("paper/paper.pdf",
                     "submission/biorxiv/pca_topology_crossdomain_replication_biorxiv.pdf"):
        path = ROOT / relative
        require(path.is_file(), f"sealed PDF absent: {relative}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        require(digest == ORIGINAL_PDF_SHA256,
                f"original sealed PDF unexpectedly changed: {relative}")

    print("ARCHIVE PROVENANCE: PASS — v0.3.0 distinct, historic DOIs preserved, original PDFs sealed")
    print("OFFLINE HASH CHECK: PASS — see the separate public Zenodo ZIP verification for external evidence")


if __name__ == "__main__":
    main()
