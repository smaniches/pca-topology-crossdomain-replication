"""Fail-closed, atomic checkpoint handling for the ablation sweep.

This module does not run statistical experiments or interpret their outputs.
The caller is responsible for the actual numerical computation.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Callable

# Required by the downstream aggregation/reporting stages. A partially written
# record must never be mistaken for an analytically usable result.
REQUIRED_RESULT_FIELDS = (
    "real_raw_max_pers", "real_pca_max_pers", "real_pca_delta", "hvg_genes"
)


def load_records(path: str) -> list[dict]:
    """Load the checkpoint, rejecting malformed or duplicate-tag records."""
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as file:
        rows = json.load(file)
    if not isinstance(rows, list):
        raise ValueError("Sweep checkpoint must contain a list of records")
    seen: set[str] = set()
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or not isinstance(row.get("tag"), str) or not row["tag"]:
            raise ValueError(f"Invalid sweep checkpoint record at index {index}")
        if row["tag"] in seen:
            raise ValueError(f"Duplicate sweep checkpoint tag: {row['tag']}")
        seen.add(row["tag"])
    return rows


def is_complete(record: dict | None, config: dict, results_dir: str) -> bool:
    """Only reuse a successful result for the identical effective configuration."""
    if not record or record.get("FAILED"):
        return False
    if any(record.get(key) != value for key, value in config.items()):
        return False
    if any(key not in record for key in REQUIRED_RESULT_FIELDS):
        return False
    artifact = Path(results_dir) / f"nulldist_{config['tag']}.pkl"
    return artifact.is_file() and artifact.stat().st_size > 0


def save_records(path: str, rows: list[dict]) -> None:
    """Atomically replace the JSON checkpoint, leaving the old file on error."""
    destination = Path(path)
    tmp_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=destination.parent,
            prefix=".sweep-checkpoint-", suffix=".tmp", delete=False,
        ) as file:
            tmp_name = file.name
            json.dump(rows, file, indent=2, default=str)
            file.flush()
            os.fsync(file.fileno())
        os.replace(tmp_name, destination)
    finally:
        if tmp_name is not None and os.path.exists(tmp_name):
            os.unlink(tmp_name)


def execute_configs(
    configs: list[dict],
    rows: list[dict],
    results_dir: str,
    checkpoint_path: str,
    run_one: Callable[[dict], dict],
    report_failure: Callable[[dict, Exception], None],
) -> tuple[list[str], int]:
    """Resume requested configs, retry failures, and persist each outcome.

    Returns (failed_tags, skipped_success_count). The caller must exit nonzero
    if failed_tags is nonempty; it must not report a completed sweep otherwise.
    """
    by_tag = {row["tag"]: index for index, row in enumerate(rows)}
    if len(by_tag) != len(rows):
        raise ValueError("Sweep checkpoint contains duplicate tags")
    requested = [config["tag"] for config in configs]
    if len(set(requested)) != len(requested):
        raise ValueError("Requested sweep configurations contain duplicate tags")

    failed: list[str] = []
    skipped = 0
    for config in configs:
        tag = config["tag"]
        index = by_tag.get(tag)
        previous = rows[index] if index is not None else None
        if is_complete(previous, config, results_dir):
            skipped += 1
            continue

        try:
            result = run_one(config)
            if not isinstance(result, dict) or result.get("tag") != tag:
                raise ValueError(f"Sweep result has an invalid tag for {tag}")
            if not is_complete(result, config, results_dir):
                raise ValueError(f"Sweep result or null-distribution artifact is incomplete: {tag}")
        except Exception as error:
            report_failure(config, error)
            result = {**config, "FAILED": True, "error": str(error)}
            failed.append(tag)

        if index is None:
            by_tag[tag] = len(rows)
            rows.append(result)
        else:
            rows[index] = result
        save_records(checkpoint_path, rows)
    return failed, skipped
