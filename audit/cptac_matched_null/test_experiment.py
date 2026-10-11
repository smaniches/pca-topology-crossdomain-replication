"""Unit and synthetic topology controls for CPTAC-TUMOR-MATCHED-NULL-002."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA

SOURCE = Path(__file__).with_name("experiment.py")
spec = importlib.util.spec_from_file_location("cptac_matched_null_for_tests", SOURCE)
experiment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(experiment)


class FeatureShuffleControls(unittest.TestCase):
    def test_exact_single_feature_empirical_distributions(self):
        rng = np.random.default_rng(831)
        x = rng.normal(size=(27, 19))
        x[:, 3] = np.arange(len(x))
        y = experiment.permute_features_independently(x, np.random.default_rng(31))
        np.testing.assert_array_equal(np.sort(y, axis=0), np.sort(x, axis=0))
        np.testing.assert_allclose(y.mean(axis=0), x.mean(axis=0), atol=1e-13)
        np.testing.assert_allclose(y.var(axis=0), x.var(axis=0), atol=1e-13)
        self.assertFalse(np.array_equal(x, y))

    def test_random_seed_determinism_and_columnwise_independence(self):
        x = np.arange(24*17,dtype=float).reshape(24,17)
        a = experiment.permute_features_independently(x,np.random.default_rng(77))
        b = experiment.permute_features_independently(x,np.random.default_rng(77))
        c = experiment.permute_features_independently(x,np.random.default_rng(78))
        np.testing.assert_array_equal(a,b)
        self.assertFalse(np.array_equal(a,c))
        self.assertFalse(np.array_equal(a[:,0].argsort(),a[:,1].argsort()))

    def test_covariance_changes_in_correlated_synthetic_data(self):
        rng=np.random.default_rng(35)
        a=rng.normal(size=(68,1))
        x=np.hstack([a+0.04*rng.normal(size=(68,1)) for _ in range(9)])
        y=experiment.permute_features_independently(x,np.random.default_rng(13))
        orig=np.corrcoef(x,rowvar=False)
        changed=np.corrcoef(y,rowvar=False)
        self.assertGreater(float(np.linalg.norm(orig-changed)),3.0)
        np.testing.assert_array_equal(np.sort(y,axis=0),np.sort(x,axis=0))

    def test_full_svd_pca_metric_parity_on_multiple_shapes(self):
        rng=np.random.default_rng(2026)
        for n,p,k in ((24,37,9),(36,500,10),(110,2000,50)):
            with self.subTest(n=n,p=p,k=k):
                x=rng.normal(size=(n,p))
                scores=experiment.exact_pca_scores(x,k)
                err=experiment.pca_distance_parity(x,scores,1e-8)
                self.assertLess(err,1e-8)
                y=experiment.permute_features_independently(x,rng)
                scores_y=experiment.exact_pca_scores(y,k)
                self.assertLess(experiment.pca_distance_parity(y,scores_y,1e-8),1e-8)
                self.assertEqual(scores_y.shape,(n,k))

    def test_invalid_matrices_rejected(self):
        for x in (np.zeros(4),np.array([[np.nan]*6]*9),np.zeros((3,5))):
            with self.assertRaises(ValueError):
                experiment.permute_features_independently(x,np.random.default_rng(2))
        with self.assertRaises(ValueError):
            experiment.exact_pca_scores(np.zeros((7,8)),12)
        with self.assertRaises(ValueError):
            experiment.exact_pca_scores(np.zeros((7,8)),5)

    def test_monte_carlo_interval_and_threshold(self):
        self.assertEqual(experiment.clopper_pearson(0,499)[0],0.0)
        self.assertEqual(experiment.clopper_pearson(499,499)[1],1.0)
        lo,hi=experiment.clopper_pearson(45,499)
        self.assertLess(lo,45/499)
        self.assertGreater(hi,45/499)
        self.assertAlmostEqual((1+45)/500,0.092)

    def test_ring_line_h1_sanity(self):
        source=experiment.haar_source()
        t=2*np.pi*np.arange(44)/44
        ring=np.column_stack((np.cos(t),np.sin(t)))
        line=np.column_stack((np.linspace(-1,1,44),np.zeros(44)))
        self.assertGreater(source.h1_max(ring),0.2)
        self.assertLess(source.h1_max(line),1e-5)

    def test_h1_unchanged_by_orthogonal_score_basis(self):
        source=experiment.haar_source()
        rng=np.random.default_rng(8)
        points=rng.normal(size=(30,6))
        q,r=np.linalg.qr(rng.normal(size=(6,6)))
        rotated=points@q
        np.testing.assert_allclose(pdist(points),pdist(rotated),rtol=1e-12)
        self.assertAlmostEqual(source.h1_max(points),source.h1_max(rotated),places=4)

    def test_existing_499_draws_are_hash_bound(self):
        prior=experiment.archived_reference()
        self.assertEqual(prior["draws"],499)
        self.assertEqual(prior["null_exceedances"],45)
        self.assertEqual(prior["p_one_sided_plus_one"],0.092)


if __name__=="__main__":
    unittest.main()
