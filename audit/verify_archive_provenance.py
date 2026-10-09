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

    require(
        bool(re.search(r'^doi:\s*["\']?10\.5281/zenodo\.21287944["\']?\s*$', citation, re.M)),
        "original DOI missing or changed in CITATION.cff",
    )
    require(ORIGINAL_DOI in readme, "original DOI omitted from README")
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

    print("ARCHIVE PROVENANCE: PASS — original DOI cited, version notice linked, PDFs sealed")
    print("EXTERNAL ZENODO RECORD: UNKNOWN — this check cannot verify archive synchronization")


if __name__ == "__main__":
    main()
