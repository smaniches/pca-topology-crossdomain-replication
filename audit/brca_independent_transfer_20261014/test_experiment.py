"""Synthetic, deliberately non-biological checks for BRCA-TRANSFER-001."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

TARGET = Path(__file__).with_name("experiment.py")
SPEC = importlib.util.spec_from_file_location("brca_transfer_test_target", TARGET)
lab = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lab)


class IndependentCohortControls(unittest.TestCase):
    def test_columnwise_shuffle_retains_exact_marginals(self):
        x = np.arange(75 * 16, dtype=float).reshape(75, 16)
        scrambled = lab.feature_permutation(x, np.random.default_rng(2026))
        np.testing.assert_array_equal(
            np.sort(x, axis=0), np.sort(scrambled, axis=0)
        )
        self.assertFalse(np.array_equal(x, scrambled))
        another = lab.feature_permutation(x, np.random.default_rng(2026))
        np.testing.assert_array_equal(scrambled, another)

    def test_haar_retains_sample_centered_covariance_and_means(self):
        x = np.random.default_rng(0).normal(size=(60, 20))
        rotated = lab.haar_rotate(x, np.random.default_rng(13))
        self.assertLess(lab.covariance_discrepancy(x, rotated), 1e-12)
        np.testing.assert_allclose(x.mean(axis=0), rotated.mean(axis=0), atol=1e-12)
        self.assertGreater(np.linalg.norm(x - rotated), 0.1)

    def test_gram_pca_metric_independent_full_svd(self):
        x = np.random.default_rng(10).normal(size=(63, 2000))
        scores = lab.pca_scores(x)
        self.assertEqual(scores.shape, (63, 50))
        self.assertLess(lab.assert_pca_parity(x, scores), 1e-7)

    def test_score_orthogonal_rotation_preserves_metric(self):
        scores = np.random.default_rng(4).normal(size=(30, 3))
        q, _ = np.linalg.qr(np.random.default_rng(2).normal(size=(3, 3)))
        np.testing.assert_allclose(pdist(scores), pdist(scores @ q), atol=1e-12)

    def test_circle_h1_is_positive_line_h1_zero(self):
        theta = np.linspace(0, 2 * np.pi, 44, endpoint=False)
        circle = np.stack([np.cos(theta), np.sin(theta)], axis=1)
        line = np.stack([np.linspace(-1, 1, 44), np.zeros(44)], axis=1)
        self.assertGreater(lab.h1_max(circle), 0.2)
        self.assertLess(lab.h1_max(line), 1e-5)

    def test_patient_matching_does_not_admit_normal_tissues(self):
        patients = [f"AO-A{i:03d}" for i in range(65)]
        columns = [patient + ".01TCGA" for patient in patients]
        df = pd.DataFrame(columns=columns + ["263d3f-I", "random"])
        clin = pd.DataFrame(
            {"Complete TCGA ID": ["TCGA-" + patient for patient in patients]}
        )
        result = lab.sample_mapping(df, clin)
        self.assertEqual(len(result), 65)
        self.assertNotIn("263d3f-I", [column for _, column in result])

    def test_duplicate_tumor_patient_rejected(self):
        patients = [f"AO-A{i:03d}" for i in range(65)]
        cols = [p + ".01TCGA" for p in patients]
        df = pd.DataFrame(columns=cols + [cols[0]])
        clin = pd.DataFrame(
            {"Complete TCGA ID": ["TCGA-" + p for p in patients]}
        )
        with self.assertRaises(ValueError):
            lab.sample_mapping(df, clin)

    def test_short_cohort_rejected_instead_of_substituting(self):
        df = pd.DataFrame(columns=["AO-A001.01TCGA"])
        clin = pd.DataFrame({"Complete TCGA ID": ["TCGA-AO-A001"]})
        with self.assertRaises(ValueError):
            lab.sample_mapping(df, clin)

    def test_malformed_pca_input_rejected(self):
        with self.assertRaises(ValueError):
            lab.pca_scores(np.zeros((10, 6)), 11)
        with self.assertRaises(ValueError):
            lab.pca_scores(np.zeros((20, 6)), 5)

    def test_reference_tail_and_monte_carlo_uncertainty(self):
        values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        result = lab.reference_report(values, observed=0.35, seed=1, name="test")
        self.assertEqual(result["exceedances"], 2)
        self.assertAlmostEqual(result["p_one_sided_plus_one"], 3 / 6)
        lo, hi = result["monte_carlo_exceedance_fraction_95pct_interval"]
        self.assertLess(lo, 2 / 5)
        self.assertGreater(hi, 2 / 5)

    def test_synthetic_full_null_pipeline_two_draws(self):
        x = np.random.default_rng(26).normal(size=(61, lab.N_FEATURES))
        old = lab.DRAWS
        try:
            lab.DRAWS = 2
            record = lab.confirm(x)
        finally:
            lab.DRAWS = old
        self.assertEqual(
            record["independent_feature_reference"]["draws"], 2
        )
        self.assertEqual(
            record["covariance_preserving_reference"]["draws"], 2
        )
        self.assertTrue(
            record["checks"]["feature_first_draw_exact_marginal_preservation"]
        )
        self.assertLess(
            record["checks"]["haar_first_draw_relative_covariance_error"], 1e-9
        )


if __name__ == "__main__":
    unittest.main()
