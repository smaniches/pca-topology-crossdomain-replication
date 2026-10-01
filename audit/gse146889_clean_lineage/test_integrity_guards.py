#!/usr/bin/env python3
"""Regression tests for clean-lineage integrity guards."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

import aggregate
import verify_reconciliation


def expect_runtime_error(fn, message: str) -> None:
    try:
        fn()
    except RuntimeError:
        return
    raise AssertionError(message)


def write_shard(path: Path, *, prepared_sha: str, worker_sha: str, base_seed: int) -> None:
    np.savez_compressed(
        path,
        task=np.asarray("mixed_gauss"),
        start=np.asarray(0, dtype=np.int64),
        count=np.asarray(2, dtype=np.int64),
        indices=np.asarray([0, 1], dtype=np.int64),
        values=np.asarray([1.0, 2.0], dtype=np.float64),
        base_seed=np.asarray(base_seed, dtype=np.int64),
        prepared_sha256=np.asarray(prepared_sha),
        worker_sha256=np.asarray(worker_sha),
    )


def main() -> None:
    original_counts = dict(aggregate.EXPECTED_COUNTS)
    aggregate.EXPECTED_COUNTS.clear()
    aggregate.EXPECTED_COUNTS["mixed_gauss"] = 2
    try:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shard = root / "mixed_gauss_00000_00002.npz"
            write_shard(
                shard,
                prepared_sha="prepared-ok",
                worker_sha="worker-ok",
                base_seed=aggregate.TASK_BASE_SEEDS["mixed_gauss"],
            )
            draws, _ = aggregate.load_shards(
                root,
                expected_prepared_sha256="prepared-ok",
                expected_worker_sha256="worker-ok",
            )
            assert np.array_equal(draws["mixed_gauss"], np.asarray([1.0, 2.0]))

            write_shard(
                shard,
                prepared_sha="wrong",
                worker_sha="worker-ok",
                base_seed=aggregate.TASK_BASE_SEEDS["mixed_gauss"],
            )
            expect_runtime_error(
                lambda: aggregate.load_shards(
                    root,
                    expected_prepared_sha256="prepared-ok",
                    expected_worker_sha256="worker-ok",
                ),
                "prepared hash mismatch was accepted",
            )

            write_shard(
                shard,
                prepared_sha="prepared-ok",
                worker_sha="wrong",
                base_seed=aggregate.TASK_BASE_SEEDS["mixed_gauss"],
            )
            expect_runtime_error(
                lambda: aggregate.load_shards(
                    root,
                    expected_prepared_sha256="prepared-ok",
                    expected_worker_sha256="worker-ok",
                ),
                "worker hash mismatch was accepted",
            )

            write_shard(
                shard,
                prepared_sha="prepared-ok",
                worker_sha="worker-ok",
                base_seed=999,
            )
            expect_runtime_error(
                lambda: aggregate.load_shards(
                    root,
                    expected_prepared_sha256="prepared-ok",
                    expected_worker_sha256="worker-ok",
                ),
                "base-seed mismatch was accepted",
            )

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"pc": 45, "pc": 0.997}\n', encoding="utf-8")
            try:
                verify_reconciliation.load_json_strict(duplicate)
            except ValueError:
                pass
            else:
                raise AssertionError("duplicate JSON key was accepted")
    finally:
        aggregate.EXPECTED_COUNTS.clear()
        aggregate.EXPECTED_COUNTS.update(original_counts)

    print(json.dumps({"integrity_guard_regressions": "PASS"}, sort_keys=True))


if __name__ == "__main__":
    main()
