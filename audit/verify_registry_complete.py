"""Fail closed if the live Git tree and tracked-file registries diverge.

Historical release ZIPs retain their own frozen checksums; this verifier
checks the *current checkout*, not any previously published Zenodo archive.
"""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_RE = re.compile(r"^\| `([^`]+)` \| (\d+) \|$")
CHECKSUM_RE = re.compile(r"^([0-9a-f]{64})  (.+)$")


def validate(root: Path) -> int:
    tracked_bytes = subprocess.run(
        ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True,
    ).stdout
    tracked = set(filter(None, tracked_bytes.decode("utf-8").split("\0")))
    manifest: dict[str, int] = {}
    for line in (root / "MANIFEST.md").read_text(encoding="utf-8").splitlines():
        match = MANIFEST_RE.fullmatch(line)
        if match:
            path = match.group(1)
            if path in manifest:
                raise ValueError(f"Duplicate manifest path: {path}")
            manifest[path] = int(match.group(2))
    checksums: dict[str, str] = {}
    for line in (root / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        match = CHECKSUM_RE.fullmatch(line)
        if not match:
            raise ValueError("Malformed checksum registry entry")
        sha, path = match.groups()
        if path in checksums:
            raise ValueError(f"Duplicate checksum path: {path}")
        checksums[path] = sha

    if set(manifest) != tracked:
        raise ValueError(f"Manifest paths differ from Git: "
                         f"missing={sorted(tracked - set(manifest))}, "
                         f"extra={sorted(set(manifest) - tracked)}")
    expected_checksums = tracked - {"checksums.sha256"}
    if set(checksums) != expected_checksums:
        raise ValueError(f"Checksum paths differ from Git: "
                         f"missing={sorted(expected_checksums - set(checksums))}, "
                         f"extra={sorted(set(checksums) - expected_checksums)}")

    for path in sorted(tracked):
        data = (root / path).read_bytes()
        if len(data) != manifest[path]:
            raise ValueError(f"Manifest byte count mismatch: {path}")
        if path in checksums and hashlib.sha256(data).hexdigest() != checksums[path]:
            raise ValueError(f"Checksum mismatch: {path}")
    return len(tracked)


def main() -> None:
    try:
        count = validate(ROOT)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"FILE REGISTRY ERROR: {error}") from error
    print(f"FILE REGISTRY PASS: {count} tracked files verified")


if __name__ == "__main__":
    main()
