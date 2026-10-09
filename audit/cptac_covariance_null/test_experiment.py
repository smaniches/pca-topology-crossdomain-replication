"""Invariant tests for CPTAC-TUMOR-COV-NULL-001.

All deterministic tests are executable without CPTAC biological data.
The optional ripser test exercises the PH statistic, not the null's validity
under unverified clinical-sample exchangeability assumptions.
"""
import importlib.util
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA

spec = importlib.util.spec_from_file_location(
    "cptac_haar_null", Path(__file__).with_name("experiment.py")
)
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


class HaarCovarianceInvariants(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(2026)
        self.x = self.rng.normal(size=(24, 37))
        self.x += self.rng.normal(size=(37,))
        self.k = 9

    def test_helmert_basis(self):
        u = ex.helmert_basis(len(self.x))
        self.assertEqual(u.shape, (24, 23))
        np.testing.assert_allclose(u.T @ u, np.eye(23), atol=1e-13)
        np.testing.assert_allclose(u.sum(axis=0), np.zeros(23), atol=1e-13)

    def test_haar_orthogonal_deterministic(self):
        q = ex.haar_orthogonal(23, np.random.default_rng(21))
        np.testing.assert_allclose(q.T @ q, np.eye(23), atol=1e-13)
        np.testing.assert_array_equal(q, ex.haar_orthogonal(23, np.random.default_rng(21)))
        self.assertFalse(np.array_equal(q, ex.haar_orthogonal(23, np.random.default_rng(22))))

    def test_full_covariance_and_mean_conservation(self):
        state = ex.score_state(self.x, self.k)
        q = ex.haar_orthogonal(23, np.random.default_rng(19))
        rotated = state["mu"] + state["u"] @ (q @ state["z"])
        np.testing.assert_allclose(rotated.mean(axis=0), self.x.mean(axis=0), atol=1e-12)
        z0 = self.x - self.x.mean(axis=0)
        zb = rotated - rotated.mean(axis=0)
        np.testing.assert_allclose(z0.T @ z0, zb.T @ zb, rtol=1e-12, atol=1e-11)
        np.testing.assert_allclose(np.linalg.svd(z0, compute_uv=False),
                                   np.linalg.svd(zb, compute_uv=False), rtol=1e-10, atol=1e-10)

    def test_exact_pca_distance_equivalence(self):
        state = ex.score_state(self.x, self.k)
        err = ex.verify_pca_parity(self.x, state, tol=1e-10)
        self.assertLess(err, 1e-10)
        q = ex.haar_orthogonal(23, np.random.default_rng(15))
        xb = state["mu"] + state["u"] @ (q @ state["z"])
        scoreb = (state["u"] @ (q @ state["left"])) * state["singular"][:self.k]
        reference = PCA(n_components=self.k, svd_solver="full").fit_transform(xb)
        np.testing.assert_allclose(pdist(scoreb), pdist(reference), rtol=1e-10, atol=1e-10)

    def test_score_generation_reproducible_with_moment_checks(self):
        state = ex.score_state(self.x, self.k)
        a = ex.rotated_scores(state, np.random.default_rng(5), check=True)
        b = ex.rotated_scores(state, np.random.default_rng(5), check=True)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(a.shape, (24, self.k))
        self.assertTrue(np.isfinite(a).all())
        self.assertFalse(np.array_equal(a, ex.rotated_scores(state, np.random.default_rng(6))))

    def test_invalid_matrices_rejected(self):
        for matrix in (np.ones(5), np.array([[float("nan")] * 8] * 7),
                       np.zeros((3, 4)), np.zeros((7, 0))):
            with self.assertRaises(ValueError):
                ex.score_state(matrix, self.k)
        with self.assertRaises(ValueError):
            ex.score_state(self.x, components=24)

    @unittest.skipUnless(importlib.util.find_spec("ripser"), "ripser not installed")
    def test_known_circle_and_flat_line(self):
        t = np.arange(40) * (2 * np.pi / 40)
        ring = np.column_stack((np.cos(t), np.sin(t)))
        line = np.column_stack((np.linspace(-1, 1, 40), np.zeros(40)))
        self.assertGreater(ex.h1_max(ring), 0.2)
        self.assertLess(ex.h1_max(line), 1e-5)

    @unittest.skipUnless(importlib.util.find_spec("ripser"), "ripser not installed")
    def test_orthogonal_feature_rotation_preserves_h1(self):
        rng = np.random.default_rng(9)
        scores = rng.normal(size=(27, 6))
        q = ex.haar_orthogonal(6, rng)
        self.assertAlmostEqual(ex.h1_max(scores), ex.h1_max(scores @ q), places=4)


if __name__ == "__main__":
    unittest.main()
