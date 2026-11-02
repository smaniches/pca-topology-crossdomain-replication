"""No-data tests for checkpoint resume, failure, and atomicity contracts."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sweep_checkpoint import execute_configs, is_complete, load_records, save_records


class SweepCheckpointTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root / "sweep_results.json"
        self.config = {
            "tag": "case-a", "imputation_rule": "zero", "n_hvg": 500,
            "n_pcs": 10, "n_perm_pipeline": 20, "n_draws_gaussian": 20,
        }
        self.errors: list[str] = []

    def result(self, config: dict) -> dict:
        return {
            **config, "real_raw_max_pers": 1.0, "real_pca_max_pers": 2.0,
            "real_pca_delta": 1.0, "hvg_genes": ["a", "b"],
        }

    def create_artifact(self, config: dict) -> None:
        (self.root / f"nulldist_{config['tag']}.pkl").write_bytes(b"fixture")

    def execute(self, configs: list[dict], rows: list[dict], runner):
        return execute_configs(
            configs, rows, str(self.root), str(self.out), runner,
            lambda c, exc: self.errors.append(f"{c['tag']}: {exc}"),
        )

    def test_reuse_requires_matching_record_and_artifact(self) -> None:
        row = self.result(self.config)
        self.assertFalse(is_complete(row, self.config, str(self.root)))
        self.create_artifact(self.config)
        self.assertTrue(is_complete(row, self.config, str(self.root)))
        self.assertFalse(is_complete({**row, "FAILED": True}, self.config, str(self.root)))
        self.assertFalse(is_complete({k: v for k, v in row.items() if k != "hvg_genes"},
                                     self.config, str(self.root)))
        changed = {**self.config, "n_perm_pipeline": 500}
        self.assertFalse(is_complete(row, changed, str(self.root)))

    def test_reuse_of_success_skips_execution(self) -> None:
        self.create_artifact(self.config)
        previous = [self.result(self.config)]
        self.assertEqual(([], 1), self.execute([self.config], previous,
                         lambda config: self.fail("must not rerun completed config")))
        self.assertFalse(self.out.exists())

    def test_failed_checkpoint_is_retried_in_place(self) -> None:
        previous = [{**self.config, "FAILED": True, "error": "old exception"}]
        calls = []
        def runner(config):
            calls.append(config["tag"])
            self.create_artifact(config)
            return self.result(config)
        self.assertEqual(([], 0), self.execute([self.config], previous, runner))
        self.assertEqual(["case-a"], calls)
        self.assertEqual(1, len(load_records(str(self.out))))
        self.assertFalse(load_records(str(self.out))[0].get("FAILED", False))

    def test_quick_count_change_forces_recompute_and_replaces_record(self) -> None:
        self.create_artifact(self.config)
        previous = [self.result(self.config)]
        full = {**self.config, "n_perm_pipeline": 500, "n_draws_gaussian": 500}
        def runner(config):
            self.create_artifact(config)
            return self.result(config)
        self.assertEqual(([], 0), self.execute([full], previous, runner))
        self.assertEqual(1, len(previous))
        self.assertEqual(500, load_records(str(self.out))[0]["n_perm_pipeline"])

    def test_failed_execution_is_persisted_and_reported(self) -> None:
        self.create_artifact(self.config)
        rows = [self.result(self.config)]
        changed = {**self.config, "n_perm_pipeline": 100}
        def runner(config):
            raise RuntimeError("synthetic failure")
        self.assertEqual((["case-a"], 0), self.execute([changed], rows, runner))
        stored = load_records(str(self.out))
        self.assertEqual(1, len(stored))
        self.assertTrue(stored[0]["FAILED"])
        self.assertEqual(100, stored[0]["n_perm_pipeline"])
        self.assertIn("synthetic failure", self.errors[0])

    def test_invalid_and_duplicate_checkpoints_rejected(self) -> None:
        for data in ({"tag": "x"}, [{"tag": "x"}, {"tag": "x"}], [{"error": "missing"}]):
            self.out.write_text(json.dumps(data), encoding="utf-8")
            with self.subTest(data=data), self.assertRaises(ValueError):
                load_records(str(self.out))
        with self.assertRaises(ValueError):
            self.execute([self.config, self.config], [], lambda _: self.result(self.config))

    def test_failed_atomic_write_preserves_prior_checkpoint(self) -> None:
        self.out.write_text('[{"tag":"old"}]', encoding="utf-8")
        original = self.out.read_bytes()
        with patch("sweep_checkpoint.os.replace", side_effect=OSError("simulated filesystem error")):
            with self.assertRaises(OSError):
                save_records(str(self.out), [self.result(self.config)])
        self.assertEqual(original, self.out.read_bytes())
        self.assertEqual([], list(self.root.glob(".sweep-checkpoint-*.tmp")))

    def test_successful_atomic_write_is_readable(self) -> None:
        save_records(str(self.out), [self.result(self.config)])
        self.assertEqual([self.result(self.config)], load_records(str(self.out)))
        self.assertEqual([], list(self.root.glob(".sweep-checkpoint-*.tmp")))


if __name__ == "__main__":
    unittest.main()
