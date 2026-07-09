"""
confound_audit_common.py -- shared statistical machinery for the confound-attribution
audit (TOPOLOGICA project), used identically across all 4 cohorts (GSE81089 pilot,
CPTAC-CCRCC, TCGA-LUAD, GSE146889).

The audit answers: "is the observed max-H1-persistence topological signal a genuine
tumor-associated feature, or an artifact of the tumor/normal class-mean-shift confound
(which separately produces near-ceiling classifier AUC in every cohort)?" via 4 controls:

1. within_class_decomposition -- does the tumor-only subset (no normal samples, so no
   class-mean-shift direction to exploit) reproduce the mixed-set max-H1 magnitude?
2. within_stratum_control -- bin samples by a continuous confound score (linear class-
   mean-shift projection, or classifier probability) into quartiles; check whether the
   topological z-score survives within each narrow-confound-range stratum.
3. residualization_control -- explicitly regress out the class-mean-shift direction
   (shift each sample by its class-specific mean minus the pooled grand mean), recompute
   PCA and max-H1, and check whether the signal and the classifier AUC collapse together.
4. block_bootstrap_ci -- resample-with-replacement from a fixed subset (e.g. tumor-only)
   and get a 95% CI on the "signal beyond null" quantity (observed - null mean), checking
   whether it excludes 0.

All 4 controls are cohort-agnostic given a preprocessed data matrix, a binary tumor/normal
mask, and null-generation parameters -- only data loading and the specific preprocessing
pipeline (HVG selection rule, imputation, log-transform) differ per cohort, which is why
the per-cohort driver scripts handle only that and delegate all statistics here.
"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict
from sklearn.metrics import roc_auc_score


# --------------------------------------------------------------------------
# Core topological statistic: max H1 (loop) persistence
# --------------------------------------------------------------------------

def max_h1_persistence(X, maxdim=1, backend="ripser"):
    """Max H1 (1-dimensional loop) persistence of a point cloud X (n_samples x n_features).

    backend='ripser' (default, used by GSE81089/TCGA-LUAD/GSE146889) computes directly
    on the point cloud via the `ripser` package.
    backend='gudhi' (used by CPTAC-CCRCC in the original audit) computes an explicit
    pairwise-Euclidean distance matrix and builds a Rips complex via `gudhi`. Both give
    the same H1 persistence diagram for the same point cloud (same Vietoris-Rips
    filtration); the backend choice in each cohort matches what was actually run so the
    cohort's committed numbers reproduce exactly.
    """
    if backend == "ripser":
        from ripser import ripser as _ripser
        res = _ripser(X, maxdim=maxdim)
        h1 = res["dgms"][1]
        finite = h1[np.isfinite(h1[:, 1])]
        if len(finite) == 0:
            return 0.0
        return float((finite[:, 1] - finite[:, 0]).max())
    elif backend == "gudhi":
        import gudhi
        from scipy.spatial.distance import pdist, squareform
        D = squareform(pdist(X, metric="euclidean"))
        rc = gudhi.RipsComplex(distance_matrix=D)
        st = rc.create_simplex_tree(max_dimension=max(maxdim, 2))
        st.persistence()
        diag = st.persistence()
        h1 = [(b, d) for (dim, (b, d)) in diag if dim == 1]
        if len(h1) == 0:
            return 0.0
        pers = [d - b for (b, d) in h1 if np.isfinite(d)]
        return max(pers) if pers else 0.0
    else:
        raise ValueError(f"unknown backend {backend!r}")


# --------------------------------------------------------------------------
# Null distributions
# --------------------------------------------------------------------------

def gaussian_null_maxh1(n_samples, n_features, n_draws, seed, n_pca=None, backend="ripser"):
    """iid standard-normal null, optionally re-projected through PCA(n_pca) to match
    the pipeline's dimensionality reduction step. Returns an array of n_draws max-H1 values."""
    rng = np.random.RandomState(seed)
    vals = np.empty(n_draws)
    for i in range(n_draws):
        X = rng.normal(size=(n_samples, n_features))
        if n_pca is not None:
            Xs = StandardScaler().fit_transform(X)
            X = PCA(n_components=min(n_pca, n_samples - 1, n_features), random_state=0).fit_transform(Xs)
        vals[i] = max_h1_persistence(X, backend=backend)
    return vals


def permutation_null_maxh1(X_real, n_draws, seed, backend="ripser"):
    """Pipeline-symmetric permutation null: independently permute each feature (gene/
    protein) column across samples, preserving each feature's marginal distribution but
    destroying sample-to-sample (row) structure. Operates directly in the space passed in
    (raw HVG space or already-PCA-projected space, matching whichever the caller wants
    a null for)."""
    rng = np.random.RandomState(seed)
    n, d = X_real.shape
    vals = np.empty(n_draws)
    for i in range(n_draws):
        perm = np.empty_like(X_real)
        for j in range(d):
            perm[:, j] = rng.permutation(X_real[:, j])
        vals[i] = max_h1_persistence(perm, backend=backend)
    return vals


def permutation_null_maxh1_pca_refit(X_allgenes_std, n_pcs, n_draws, seed, backend="ripser"):
    """Pipeline-symmetric permutation null that also re-fits PCA on each permuted draw
    (used for the residualization control's PCA-space null, and by GSE146889/CPTAC for
    their raw-gene permutation nulls that reselect HVG post-permutation)."""
    rng = np.random.RandomState(seed)
    n, d = X_allgenes_std.shape
    vals = np.empty(n_draws)
    for i in range(n_draws):
        perm = np.empty_like(X_allgenes_std)
        for j in range(d):
            perm[:, j] = rng.permutation(X_allgenes_std[:, j])
        pca = PCA(n_components=min(n_pcs, n - 1, d), random_state=42)
        Xp_pca = pca.fit_transform(perm)
        vals[i] = max_h1_persistence(Xp_pca, backend=backend)
    return vals


def zscore_pvalue(observed, null_draws):
    """z-score and one-sided (null >= observed) empirical p-value against a null distribution."""
    null_draws = np.asarray(null_draws)
    mean = null_draws.mean()
    sd = null_draws.std(ddof=1) if len(null_draws) > 1 else null_draws.std()
    z = (observed - mean) / sd if sd > 0 else np.nan
    p = (np.sum(null_draws >= observed) + 1) / (len(null_draws) + 1)
    return dict(null_mean=float(mean), null_sd=float(sd), z=float(z), p=float(p))


# --------------------------------------------------------------------------
# Control 1: within-class decomposition
# --------------------------------------------------------------------------

def within_class_decomposition(X_mixed, tumor_mask, n_gauss=500, n_perm=2000,
                                n_features=None, n_pca=None, backend="ripser",
                                seed_base=100):
    """Compute observed max-H1 for mixed / tumor-only / normal-only subsets, plus matched
    Gaussian and permutation nulls for each subset (matched in sample count, so that a
    smaller subset is compared against a null with the SAME n, not the full-cohort null).

    X_mixed: (n_samples, n_features) matrix already at the analysis space (raw HVG or
        PCA-projected -- whichever the cohort's `real_pca_delta`/`real_max_h1` refers to).
    tumor_mask: boolean array, True = tumor.
    n_features: number of raw features BEFORE PCA, needed for the Gaussian null shape when
        n_pca is set (i.e. X_mixed is already PCA-projected but the null should be drawn in
        raw gene space and then projected, matching the pipeline). If X_mixed IS the raw
        space (no PCA), leave n_pca=None and n_features=X_mixed.shape[1].
    """
    if n_features is None:
        n_features = X_mixed.shape[1]
    normal_mask = ~tumor_mask
    subsets = {"mixed": np.ones(len(tumor_mask), dtype=bool), "tumor": tumor_mask, "normal": normal_mask}
    out = {}
    for i, (name, mask) in enumerate(subsets.items()):
        X_sub = X_mixed[mask]
        n_sub = X_sub.shape[0]
        obs = max_h1_persistence(X_sub, backend=backend)
        gauss = gaussian_null_maxh1(n_sub, n_features, n_gauss, seed=seed_base + i * 100,
                                     n_pca=n_pca, backend=backend)
        perm = permutation_null_maxh1(X_sub, n_perm, seed=seed_base + i * 100 + 50, backend=backend)
        out[name] = dict(n=n_sub, observed=obs,
                          gauss=zscore_pvalue(obs, gauss), perm=zscore_pvalue(obs, perm),
                          gauss_draws=gauss, perm_draws=perm)
    return out


# --------------------------------------------------------------------------
# Control 2: within-stratum (quartile) control
# --------------------------------------------------------------------------

def confound_direction(X_std, tumor_mask):
    """Linear class-mean-shift direction (unit vector) in a standardized feature space:
    mean(tumor) - mean(normal), normalized. Its projection is the "confound score" used
    to bin samples for the within-stratum control."""
    mean_tumor = X_std[tumor_mask].mean(axis=0)
    mean_normal = X_std[~tumor_mask].mean(axis=0)
    diff = mean_tumor - mean_normal
    u = diff / np.linalg.norm(diff)
    return u, X_std @ u


def within_stratum_control(X_by_stratum, n_gauss=300, backend="ripser", seed_base=1000,
                            perm_source_by_stratum=None, n_perm=300):
    """For each stratum (quartile bin, given as a list of (n_samples, X_stratum) or just
    X_stratum matrices), compute observed max-H1 and a matched Gaussian null (and,
    optionally, a matched permutation null re-selected from `perm_source_by_stratum`,
    the raw pre-HVG-selection gene matrix restricted to that stratum's samples, for
    cohorts that reselect HVG genes within the permuted stratum)."""
    rows = []
    n_features = X_by_stratum[0].shape[1]
    for q, X_q in enumerate(X_by_stratum):
        n_q = X_q.shape[0]
        obs_q = max_h1_persistence(X_q, backend=backend)
        gauss_q = gaussian_null_maxh1(n_q, n_features, n_gauss, seed=seed_base + q, backend=backend)
        row = dict(quartile=f"Q{q+1}", n=n_q, observed=obs_q, **{f"gauss_{k}": v for k, v in zscore_pvalue(obs_q, gauss_q).items()})
        if perm_source_by_stratum is not None:
            perm_q = permutation_null_maxh1_pca_refit(perm_source_by_stratum[q], n_pcs=min(50, n_q - 1),
                                                       n_draws=n_perm, seed=seed_base + 5000 + q, backend=backend)
            row.update({f"perm_{k}": v for k, v in zscore_pvalue(obs_q, perm_q).items()})
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Control 3: residualization
# --------------------------------------------------------------------------

def residualize_class_mean(X_std, tumor_mask):
    """Remove the class-mean-shift confound: shift each sample by (its class mean - the
    pooled grand mean), so both classes end up with the identical (grand) mean afterward.
    This is a linear projection-based residualization, not a full linear-model regression
    -- it removes exactly the 1-D direction used to define the confound score."""
    grand_mean = X_std.mean(axis=0)
    X_resid = X_std.copy()
    mean_tumor = X_std[tumor_mask].mean(axis=0)
    mean_normal = X_std[~tumor_mask].mean(axis=0)
    X_resid[tumor_mask] = X_std[tumor_mask] - mean_tumor + grand_mean
    X_resid[~tumor_mask] = X_std[~tumor_mask] - mean_normal + grand_mean
    max_abs_diff = np.max(np.abs(X_resid[tumor_mask].mean(axis=0) - X_resid[~tumor_mask].mean(axis=0)))
    return X_resid, float(max_abs_diff)


def cv_auc(X, y, seed=42, n_splits=5):
    clf = LogisticRegression(max_iter=5000, C=1.0)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    proba = cross_val_predict(clf, X, y, cv=skf, method="predict_proba")[:, 1]
    return float(roc_auc_score(y, proba))


def residualization_control(X_std, tumor_mask, n_pcs=50, n_gauss=500, n_perm=500,
                             backend="ripser", seed_base=9000):
    """Full residualization control: (a) residualize, (b) refit PCA on the residualized
    space, (c) recompute max-H1 there, (d) recompute matched Gaussian/permutation nulls in
    PCA(n_pcs) space, (e) check CV AUC collapses toward 0.5 in the residualized space vs
    the intact (pre-residualization) space."""
    y = tumor_mask.astype(int)
    pca_intact = PCA(n_components=n_pcs, random_state=42)
    X_pca_intact = pca_intact.fit_transform(X_std)
    obs_intact = max_h1_persistence(X_pca_intact, backend=backend)
    auc_intact = cv_auc(X_pca_intact, y)

    X_resid, max_abs_diff = residualize_class_mean(X_std, tumor_mask)
    pca_resid = PCA(n_components=n_pcs, random_state=42)
    X_pca_resid = pca_resid.fit_transform(X_resid)
    obs_resid = max_h1_persistence(X_pca_resid, backend=backend)
    auc_resid = cv_auc(X_pca_resid, y)

    n = X_std.shape[0]
    n_features = X_std.shape[1]
    gauss_pca = gaussian_null_maxh1(n, n_features, n_gauss, seed=seed_base, n_pca=n_pcs, backend=backend)
    perm_pca = permutation_null_maxh1_pca_refit(X_std, n_pcs, n_perm, seed=seed_base + 1, backend=backend)

    return dict(
        max_abs_classmean_diff_after_resid=max_abs_diff,
        intact=dict(observed=obs_intact, auc=auc_intact, gauss=zscore_pvalue(obs_intact, gauss_pca),
                    perm=zscore_pvalue(obs_intact, perm_pca)),
        residualized=dict(observed=obs_resid, auc=auc_resid, gauss=zscore_pvalue(obs_resid, gauss_pca),
                           perm=zscore_pvalue(obs_resid, perm_pca)),
        gauss_draws=gauss_pca, perm_draws=perm_pca,
    )


# --------------------------------------------------------------------------
# Control 4: block-bootstrap CI
# --------------------------------------------------------------------------

def block_bootstrap_ci(X_fixed, null_mean, null_sd, n_boot=2000, seed=42, backend="ripser"):
    """Resample rows of X_fixed with replacement n_boot times, recompute max-H1 each time,
    and report the empirical 95% CI of (observed - null_mean) ["delta"] and of
    (observed - null_mean)/null_sd ["z"]. If the delta CI excludes 0, the signal-beyond-
    null is not an artifact of any single influential sample."""
    rng = np.random.default_rng(seed)
    n = X_fixed.shape[0]
    deltas = np.empty(n_boot)
    zs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb = X_fixed[idx]
        v = max_h1_persistence(Xb, backend=backend)
        deltas[b] = v - null_mean
        zs[b] = (v - null_mean) / null_sd if null_sd > 0 else np.nan
    ci_delta = (float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5)))
    ci_z = (float(np.percentile(zs, 2.5)), float(np.percentile(zs, 97.5)))
    return dict(deltas=deltas, zs=zs, ci_delta=ci_delta, ci_z=ci_z, excludes_zero=ci_delta[0] > 0)
