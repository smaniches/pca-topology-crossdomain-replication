"""Read-only byte-level audit of prior v0.2.0 and published v0.3.0 Zenodo ZIPs.

Public API and archive bytes only: no accounts, secrets, deposits, edits,
publication actions, or new release tags. Fails closed if files, hashes,
concept/version linkage or the v0.3.0 registered source files disagree.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import urllib.parse
import urllib.request
import zipfile

OLD_RECORD = "23272317"
NEW_RECORD = None  # discovered from public probe evidence; never guessed
EXPECTED_CONCEPT = "21287943"
EXPECTED_VERSION = "v0.3.0"
SOURCE_COMMIT = "bb453a01874a97d996bdd19a7eb571c3c2bc82ae"
PDF_SHA256 = "cfc2410130c5872979a62e08cecfb964aa5412d42192c0c632c22d48e1983aac"
MAX_ZIP_BYTES = 30 * 1024 * 1024
USER_AGENT = "pca-topology-zenodo-archive-verifier/1.0"


def fetch(url: str, max_bytes: int) -> bytes:
    if not url.startswith("https://zenodo.org/"):
        raise ValueError(f"Unexpected initial archive source: {url}")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=35) as res:
        body = res.read(max_bytes + 1)
    if len(body) > max_bytes:
        raise ValueError("Archive exceeds declared safe download bound")
    return body


def record(record_id: str) -> dict:
    data = json.loads(fetch(f"https://zenodo.org/api/records/{record_id}", 3_000_000))
    if str(data.get("id")) != record_id:
        raise ValueError("Zenodo API returned wrong record")
    if str(data.get("conceptrecid")) != EXPECTED_CONCEPT:
        raise ValueError("Archive is in unexpected DOI version family")
    return data


def get_zip_file(r: dict) -> tuple[dict, bytes]:
    files = r.get("files")
    if not isinstance(files, list):
        raise ValueError("Zenodo API did not list downloadable file metadata")
    archive = next((f for f in files
                    if str(f.get("key", "")).endswith(".zip")), None)
    if archive is None:
        raise ValueError("Zip archive missing from Zenodo record")
    if archive.get("size", 0) > MAX_ZIP_BYTES:
        raise ValueError("Archive is larger than safe bound")
    links = archive.get("links") or {}
    candidates = [
        links.get("self"),
        links.get("download"),
        f"https://zenodo.org/api/records/{r['id']}/files/"
        + urllib.parse.quote(str(archive["key"]), safe="") + "/content",
    ]
    last = None
    for link in candidates:
        if not link:
            continue
        try:
            body = fetch(link, MAX_ZIP_BYTES)
            declared_size = int(archive["size"])
            if len(body) != declared_size:
                raise ValueError(f"ZIP size differs: {len(body)} != {declared_size}")
            checksum = archive.get("checksum", "")
            if not checksum.startswith("md5:"):
                raise ValueError(f"Unsupported file checksum format: {checksum}")
            actual_md5 = hashlib.md5(body).hexdigest()
            if actual_md5 != checksum.split(":", 1)[1]:
                raise ValueError("Zenodo archive bytes do not match its MD5 checksum")
            return archive, body
        except Exception as exc:
            last = str(exc)
    raise ValueError("Could not verify Zenodo ZIP: " + str(last))


def members(z: zipfile.ZipFile) -> dict[str, bytes]:
    result = {}
    for info in z.infolist():
        if info.is_dir():
            continue
        # A GitHub-generated source ZIP contains one enclosing project
        # directory, including a commit prefix. Do not trust arbitrary paths.
        bits = Path(info.filename).parts
        if len(bits) < 2 or any(bit == ".." for bit in bits):
            raise ValueError("Unexpected path structure inside archival ZIP")
        rel = "/".join(bits[1:])
        if rel in result:
            raise ValueError("Duplicate archive path: " + rel)
        result[rel] = z.read(info)
    return result


def verify_sealed_pdf(m: dict, label: str) -> None:
    paths = [
        "paper/paper.pdf",
        "submission/biorxiv/pca_topology_crossdomain_replication_biorxiv.pdf",
    ]
    for p in paths:
        if p not in m:
            raise ValueError(f"{label}: original PDF missing: {p}")
        if hashlib.sha256(m[p]).hexdigest() != PDF_SHA256:
            raise ValueError(f"{label}: original PDF digest changed: {p}")


def inspect_old(r: dict, raw: bytes) -> dict:
    if str(r.get("metadata", {}).get("version")) != "v0.2.0":
        raise ValueError("Previous Zenodo record does not identify v0.2.0")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        m = members(z)
    verify_sealed_pdf(m, "Previous v0.2.0")
    return {"version": "v0.2.0", "files": len(m),
            "original_pdf_verified": True}


def inspect_new(r: dict, raw: bytes) -> dict:
    if str(r.get("metadata", {}).get("version")) != EXPECTED_VERSION:
        raise ValueError("New Zenodo record has an unexpected version")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        m = members(z)
    verify_sealed_pdf(m, "October v0.3.0")
    required = [
        "CITATION.cff",
        "release/v0.3.0/RELEASE_NOTES.md",
        "paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md",
        "audit/cptac_covariance_null/PROTOCOL.md",
        "results/cptac_resid_null_sensitivity_20261009/499_draws.json",
        "results/cptac_tumor_covariance_null_20261009/499_draws.json",
        "results/cptac_tumor_covariance_null_20261009/REPORT.md",
        "audit/cptac_matched_null/PROTOCOL.md",
        "audit/cptac_matched_null/verify_evidence.py",
        "results/cptac_tumor_matched_null_20261010/REPORT.md",
        "results/cptac_tumor_matched_null_20261010/499_draws.json",
        "docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md",
        "checksums.sha256",
        "MANIFEST.md",
    ]
    absent = [p for p in required if p not in m]
    if absent:
        raise ValueError("v0.3.0 required release contents missing: " + str(absent))
    citation = m["CITATION.cff"].decode("utf-8")
    if "version: 0.3.0" not in citation:
        raise ValueError("Archive does not carry v0.3.0 citation metadata")
    registry = m["checksums.sha256"].decode("utf-8").splitlines()
    checks = {}
    for line in registry:
        if not line.strip():
            continue
        if len(line) < 67 or line[64:66] != "  ":
            raise ValueError("Malformed archived checksum registry entry")
        digest, path = line[:64], line[66:]
        if path not in m:
            raise ValueError("Hash-registered file missing: " + path)
        actual = hashlib.sha256(m[path]).hexdigest()
        if actual != digest:
            raise ValueError("Archived file failed its registered SHA256: " + path)
        checks[path] = digest
    if not checks or len(checks) + 1 != len(m):
        raise ValueError("Registered inventory differs from exact ZIP file set")
    matched = "results/cptac_tumor_matched_null_20261010/499_draws.json"
    if (hashlib.sha256(m[matched]).hexdigest() !=
            "7f64f411c840561fb9468a7edf887836b4a9b93b5d838ad386f43211f5e87b09"):
        raise ValueError("Matched-null raw 499-draw evidence changed in the archive")
    if ("docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md" not in m or
            b"NOT STARTED" not in m["docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md"]):
        raise ValueError("New independent-cohort phase was not kept separate")

    return {"version": EXPECTED_VERSION, "files": len(m),
            "registered_files_verified": len(checks),
            "original_pdf_verified": True,
            "required_files_present": required,
            "release_notes_bytes": len(m["release/v0.3.0/RELEASE_NOTES.md"])}


def main():
    out = Path(os.environ.get("ZENODO_BYTES_REPORT",
                              "artifacts/zenodo_v0.2.0_file_integrity.json"))
    report = {
        "status": "UNKNOWN", "source_commit": SOURCE_COMMIT,
        "original_doi": "10.5281/zenodo.21287944",
        "previous_doi": "10.5281/zenodo.23272317",
        "new_doi": None,
        "concept_doi": "10.5281/zenodo.21287943",
        "archive_records": [],
    }
    try:
        probe_path = Path(os.environ.get(
            "ZENODO_PROBE_FILE", "artifacts/zenodo_v0.3.0_probe.json"
        ))
        probe = json.loads(probe_path.read_text(encoding="utf-8"))
        if (probe.get("status") != "VERIFIED" or
                probe.get("family_match_with_old") is not True):
            raise ValueError("v0.3.0 DOI and version family have not been verified")
        doi = probe.get("version_doi") or ""
        if not doi.startswith("10.5281/zenodo."):
            raise ValueError("Unexpected Zenodo DOI identifier")
        new_record = doi.rsplit(".", 1)[1]
        if not new_record.isdigit() or new_record == OLD_RECORD:
            raise ValueError("Invalid or reused v0.3.0 version DOI")
        report["new_doi"] = doi
        for rid, inspector in [(OLD_RECORD, inspect_old),
                               (new_record, inspect_new)]:
            r = record(rid)
            meta, contents = get_zip_file(r)
            facts = inspector(r, contents)
            facts.update({
                "record_id": rid,
                "doi": r["doi"],
                "conceptrecid": str(r["conceptrecid"]),
                "zenodo_file": meta["key"],
                "zenodo_declared_md5": meta["checksum"],
                "archive_zip_sha256": hashlib.sha256(contents).hexdigest(),
                "archive_bytes": len(contents),
            })
            report["archive_records"].append(facts)
        report["status"] = "VERIFIED"
    except Exception as exc:
        report["status"] = "UNKNOWN_OR_FAILED"
        report["error"] = str(exc)[:600]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("ZENODO ARCHIVE BYTE VERIFICATION:", report["status"])
    for archive in report["archive_records"]:
        print("  Archive:", archive["doi"], archive["zenodo_file"],
              archive["archive_bytes"], "bytes")
    if report.get("error"):
        print("  Blocker:", report["error"])
    print("  Machine-readable evidence:", out)


if __name__ == "__main__":
    main()
