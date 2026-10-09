"""Non-regression tests for the new CPTAC sensitivity analysis only."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("experiment.py")
spec = importlib.util.spec_from_file_location("cptac_null_experiment", MODULE_PATH)
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


class ConditionalPermutationTests(unittest.TestCase):
    def setUp(self):
        self.mask = np.array([True] * 8 + [False] * 8)
        rng = np.random.default_rng(7)
        self.X = rng.normal(size=(16, 5))
        self.X[self.mask, :] += np.array([3, 2, 0, 0, 1])
        self.X[~self.mask, :] -= np.array([2, 0, 1, 0, 1])

    def test_preserves_all_conditional_column_values(self):
        sur = ex.shuffle_columns_within_classes(self.X, self.mask, np.random.default_rng(17))
        for cls in (True, False):
            for col in range(self.X.shape[1]):
                np.testing.assert_array_equal(np.sort(sur[self.mask == cls, col]),
                                              np.sort(self.X[self.mask == cls, col]))
        self.assertFalse(np.array_equal(self.X, sur))

    def test_reproducible_same_seed_and_different_seed(self):
        a = ex.shuffle_columns_within_classes(self.X, self.mask, np.random.default_rng(1))
        b = ex.shuffle_columns_within_classes(self.X, self.mask, np.random.default_rng(1))
        c = ex.shuffle_columns_within_classes(self.X, self.mask, np.random.default_rng(2))
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(a, c))

    def test_class_means_zero_after_transform_for_real_and_null(self):
        sur = ex.shuffle_columns_within_classes(self.X, self.mask, np.random.default_rng(55))
        for matrix in (self.X, sur):
            residual, diff = ex.residualize_class_mean(matrix, self.mask)
            self.assertLess(diff, 1e-12)
            np.testing.assert_allclose(residual[self.mask].mean(axis=0),
                                       residual[~self.mask].mean(axis=0), atol=1e-12)
            pcs = ex.residualized_pca(matrix, self.mask, n_pcs=3)
            self.assertEqual(pcs.shape, (16, 3))
            self.assertTrue(np.isfinite(pcs).all())

    def test_reject_invalid_inputs(self):
        with self.assertRaises(ValueError):
            ex.shuffle_columns_within_classes(self.X, np.ones(16, dtype=bool),
                                              np.random.default_rng(5))
        bad = self.X.copy()
        bad[0, 0] = np.nan
        with self.assertRaises(ValueError):
            ex.shuffle_columns_within_classes(bad, self.mask, np.random.default_rng(5))

    @unittest.skipUnless(importlib.util.find_spec("ripser"), "ripser not installed")
    def test_synthetic_ring_vs_line_and_class_shift_invariance(self):
        n = 36
        t = 2 * np.pi * np.arange(n) / n
        ring = np.column_stack((np.cos(t), np.sin(t)))
        line = np.column_stack((np.linspace(-1, 1, n), np.zeros(n)))
        mask = np.arange(n) % 2 == 0
        ring_with_class_shift = ring.copy()
        ring_with_class_shift[mask, 0] += 5
        ring_with_class_shift[~mask, 0] -= 5
        ring_h1 = ex.statistic(ring, mask)
        line_h1 = ex.statistic(line, mask)
        shift_h1 = ex.statistic(ring_with_class_shift, mask)
        self.assertGreater(ring_h1, 0.2)
        self.assertLess(line_h1, 1e-6)
        self.assertAlmostEqual(ring_h1, shift_h1, places=5)

    @unittest.skipUnless(importlib.util.find_spec("ripser"), "ripser not installed")
    def test_small_deterministic_statistic(self):
        rng = np.random.default_rng(123)
        X = rng.normal(size=(22, 5))
        mask = np.arange(22) % 2 == 0
        self.assertEqual(ex.statistic(X, mask), ex.statistic(X, mask))


if __name__ == "__main__":
    unittest.main()
