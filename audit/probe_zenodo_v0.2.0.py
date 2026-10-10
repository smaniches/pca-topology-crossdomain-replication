"""Read-only reconciliation of Zenodo v0.2.0 against a published GitHub release.

This script never creates a deposit, edits metadata, or assumes that GitHub
Release publication proves Zenodo ingestion. It probes public Zenodo APIs
and the repository DOI badge from an internet-connected GitHub runner.
Results are written as JSON even if unavailable or still pending.

Only status VERIFIED means a publicly retrievable v0.2.0 record with a
matching GitHub repository association and a published file was observed.
"""
import datetime
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO = "smaniches/pca-topology-crossdomain-replication"
GITHUB_REPO_ID = 1293670390
OLD_ZENODO_RECORD_ID = "21287944"
GITHUB_RELEASE = f"https://github.com/{REPO}/releases/tag/v0.2.0"
USER_AGENT = "research-archive-provenance-check/1.0 (public read-only)"


def request(url, timeout=14):
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            raw = response.read(3_000_000)
            body = None
            try:
                body = json.loads(raw)
            except (ValueError, UnicodeDecodeError):
                pass
            return {
                "ok": True, "status": response.status,
                "final_url": response.url, "json": body,
            }
    except Exception as exc:
        code = exc.code if isinstance(exc, urllib.error.HTTPError) else None
        return {"ok": False, "status": code, "error": str(exc)[:300]}


def records_from(data):
    if isinstance(data, dict) and isinstance(data.get("hits"), dict):
        return data["hits"].get("hits", [])
    if isinstance(data, dict) and data.get("id") and data.get("metadata"):
        return [data]
    return []


def file_list(record):
    files = record.get("files")
    if isinstance(files, list):
        return [{"name": f.get("key") or f.get("filename"),
                 "checksum": f.get("checksum"), "size": f.get("size")}
                for f in files]
    if isinstance(files, dict):
        return [{"name": k, "checksum": v.get("checksum"), "size": v.get("size")}
                for k, v in files.get("entries", {}).items()]
    return []


def summarize_record(record):
    meta = record.get("metadata", {})
    version = str(meta.get("version") or "")
    linked = json.dumps(
        {
            "related_identifiers": meta.get("related_identifiers"),
            "identifiers": meta.get("identifiers"),
            "links": record.get("links"),
            "description": meta.get("description"),
            "notes": meta.get("notes"),
        }, sort_keys=True, default=str,
    ).lower()
    return {
        "id": str(record.get("id") or ""),
        "doi": record.get("doi") or meta.get("doi"),
        "conceptdoi": record.get("conceptdoi"),
        "conceptrecid": str(record.get("conceptrecid") or ""),
        "title": meta.get("title"),
        "version": version,
        "github_link_confirmed": REPO.lower() in linked,
        "version_matches": version.lstrip("v") == "0.2.0",
        "files": file_list(record),
        "record_url": (record.get("links") or {}).get("html"),
    }


def main():
    out = Path(os.getenv("PROVENANCE_OUTPUT", "artifacts/zenodo_v0.2.0_probe.json"))
    evidence = {
        "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_release": GITHUB_RELEASE,
        "historical_doi": "10.5281/zenodo.21287944",
        "status": "UNKNOWN",
        "attempts": [],
        "matching_records": [],
        "family_match_with_old": None,
    }
    sources = [
        ("historical", f"https://zenodo.org/api/records/{OLD_ZENODO_RECORD_ID}"),
        ("badge", f"https://zenodo.org/badge/latestdoi/{GITHUB_REPO_ID}"),
        ("repository_query", "https://zenodo.org/api/records?q=" +
         urllib.parse.quote(REPO) + "&size=40"),
        ("title_query", "https://zenodo.org/api/records?q=" +
         urllib.parse.quote('"Does PCA Reveal or Fabricate Topological Structure"') +
         "&size=40"),
    ]
    found = {}
    old_concept = None
    reachable = False
    for round_idx in range(3):
        for kind, url in sources:
            res = request(url)
            reachable = reachable or res["ok"]
            evidence["attempts"].append({
                "round": round_idx + 1, "source": kind,
                "status": res.get("status"), "ok": res["ok"],
                "url": url, "final_url": res.get("final_url"),
                "error": res.get("error"),
            })
            if kind == "badge" and res["ok"]:
                match = re.search(r"(?:zenodo\.org/(?:records/)?|zenodo\.)(\d+)", res["final_url"])
                if match:
                    extra = request("https://zenodo.org/api/records/" + match.group(1))
                    if extra.get("ok"):
                        for record in records_from(extra.get("json")):
                            found[str(record.get("id"))] = record
            for record in records_from(res.get("json")):
                found[str(record.get("id"))] = record
                if kind == "historical":
                    old_concept = str(record.get("conceptrecid") or record.get("id"))
        candidates = [summarize_record(v) for v in found.values()]
        matches = [c for c in candidates
                   if c["version_matches"] and c["github_link_confirmed"]]
        if any(c["files"] and c["doi"] for c in matches):
            break
        if round_idx < 2:
            time.sleep(8)

    evidence["matching_records"] = [
        c for c in [summarize_record(v) for v in found.values()]
        if c["github_link_confirmed"] or c["version_matches"]
    ]
    verified = [c for c in evidence["matching_records"]
                if c["version_matches"] and c["github_link_confirmed"]
                and c["doi"] and c["files"]]
    if verified:
        evidence["status"] = "VERIFIED"
        evidence["version_doi"] = verified[0]["doi"]
        evidence["record_url"] = verified[0]["record_url"]
        if old_concept:
            evidence["family_match_with_old"] = (
                verified[0]["conceptrecid"] == old_concept
            )
    elif reachable:
        evidence["status"] = "PENDING_OR_NOT_FOUND"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print("ZENODO VERIFICATION:", evidence["status"])
    print("Historical version-family relationship:",
          evidence["family_match_with_old"])
    if evidence.get("version_doi"):
        print("New DOI:", evidence["version_doi"])
    print("Machine-readable evidence:", out)


if __name__ == "__main__":
    main()
