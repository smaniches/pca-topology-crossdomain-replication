#!/usr/bin/env python3
"""
Reproduces the GSE146889 confirmatory replication cohort end-to-end from the
public GEO accession (no local checkpoints required): fetch raw data,
preprocess, compute persistent homology (real data, both raw-HVG and PCA(50)
spaces), run both null models (Gaussian and pipeline-symmetric permutation),
compute the intrinsic-dimension gate and the tumor/normal CV-AUC confound
check, and write final_results_table_GSE146889.csv.

Cohort: GSE146889, bulk RNA-seq (gene-level RPKM), mismatch-repair-deficiency
study, 176 samples (91 tumor, 85 paired-normal; colorectal/endometrial/
ovarian). This is the 2nd pre-registered confirmatory test in the TOPOLOGICA
PCA-topology program (see ../../prereg/PREREGISTRATION.md for the locked
protocol; the GSE81089 NSCLC pilot is excluded from pre-registration).

Expected results (tolerance 1e-3 unless noted; see results/replication_GSE146889/
report_GSE146889.md Section 5 for the full table):
    real max-H1 persistence, raw HVG (top-2000 by variance)  = 6.204
    real max-H1 persistence, PCA(50)                          = 7.564
    5-fold CV AUC (tumor vs normal, PCA50 features)            = 0.925 (+/- 0.02)
    intrinsic dimension (MLE estimator), PCA(50) space         = 8.53

Usage:
    python3 01_fetch_preprocess_and_results_table.py
        # full run: N_GAUSS=500, N_PERM=2000 (as pre-registered).
        # On a single modern core this takes on the order of 30-60 minutes,
        # dominated by the 2000 permutation-null persistent-homology computations.

    python3 01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
        # fast smoke test (a few minutes) -- NOT the pre-registered draw count,
        # for verifying the pipeline runs end-to-end only. Null distributions
        # from n=20 draws are noisy; z-scores/p-values from a smoke test should
        # not be treated as the confirmatory result.

Data source: this script re-fetches the raw count/RPKM matrix directly from
NCBI GEO's public FTP mirror, the same source used in the original replication
analysis (documented explicitly here, not assumed):
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/suppl/GSE146889_GeneCount.tsv.gz
"""
import argparse
import os
import time
import warnings
import urllib.request

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from ripser import ripser

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HERE = os.path.dirname(os.path.abspath(__file__))

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/suppl/GSE146889_GeneCount.tsv.gz"
GEO_LOCAL = os.path.join(_HERE, "GSE146889_GeneCount.tsv.gz")

# Expected values from the confirmatory-test report (results/replication_GSE146889/report_GSE146889.md)
EXPECTED = {
    "real_max_h1_raw_hvg": 6.204,
    "real_max_h1_pca50": 7.564,
    "cv_auc_mean": 0.925,
    "d_int_mle_pca50": 8.53,
}
TOL = 1e-3          # tolerance for the persistence-diagram numbers (deterministic given fixed seeds)
TOL_CV_AUC = 0.02   # CV-AUC has fold-assignment variance; tolerance per task spec
TOL_DIM = 0.1        # intrinsic-dim estimate tolerance


def fetch_data():
    if not os.path.isfile(GEO_LOCAL):
        print(f"Downloading {GEO_URL} ...")
        urllib.request.urlretrieve(GEO_URL, GEO_LOCAL)
    else:
        print(f"Using cached download at {GEO_LOCAL}")
    print("File size (bytes):", os.path.getsize(GEO_LOCAL))


def load_and_preprocess():
    df_full = pd.read_csv(GEO_LOCAL, sep="\t")

    rpkm_cols = [c for c in df_full.columns if c.endswith("_rpkm")]
    meta_cols = ["GeneId", "GeneName", "GeneBiotype"]
    rpkm_df = df_full[meta_cols + rpkm_cols].copy()

    # Drop genes with any NaN across samples (data quirk: a small number of
    # gene rows carry NaN RPKM values; documented, not silently imputed).
    nan_row_mask = rpkm_df[rpkm_cols].isna().any(axis=1)
    print(f"Dropping {nan_row_mask.sum()} gene rows with NaN RPKM values (of {len(rpkm_df)})")
    rpkm_clean = rpkm_df[~nan_row_mask].copy()

    X_raw = rpkm_clean[rpkm_cols].T.copy()  # samples x genes
    X_raw.columns = rpkm_clean["GeneId"].values
    sample_names = [c.replace("_rpkm", "") for c in rpkm_cols]
    X_raw.index = sample_names
    labels = np.array(["tumor" if "_tumor_" in s else "normal" for s in sample_names])
    print("Sample matrix:", X_raw.shape, "| tumor:", (labels == "tumor").sum(), "normal:", (labels == "normal").sum())

    np.random.seed(42)
    X_log = np.log1p(X_raw.values.astype(np.float64))

    gene_var = X_log.var(axis=0)
    nonzero_var_mask = gene_var > 0
    print("Genes retained (nonzero variance):", nonzero_var_mask.sum(), "of", len(gene_var))
    X_nz = X_log[:, nonzero_var_mask]
    gene_var_nz = gene_var[nonzero_var_mask]

    order = np.argsort(gene_var_nz)  # ascending
    bottom_idx = order[:2000]
    top_idx = order[-2000:]

    X_hvg = X_nz[:, top_idx]
    X_bottom = X_nz[:, bottom_idx]
    print("HVG top-2000 shape:", X_hvg.shape, "| bottom-2000 shape:", X_bottom.shape)

    return X_nz, X_hvg, X_bottom, labels


def standardize(X):
    return StandardScaler().fit_transform(X)


def max_h1_persistence(dgms):
    h1 = dgms[1]
    if len(h1) == 0:
        return 0.0
    finite = h1[np.isfinite(h1[:, 1])]
    if len(finite) == 0:
        return 0.0
    pers = finite[:, 1] - finite[:, 0]
    return float(pers.max())


def compute_ph(X, maxdim=2):
    return ripser(X, maxdim=maxdim)["dgms"]


def pipeline_max_h1(X_gene_space, n_pca=50, seed=42):
    """Standardize -> compute raw-space max-H1 -> PCA(n_pca) -> compute PCA-space max-H1."""
    Xs = StandardScaler().fit_transform(X_gene_space)
    dgms_raw = ripser(Xs, maxdim=1)["dgms"]
    raw_stat = max_h1_persistence(dgms_raw)
    pca = PCA(n_components=n_pca, random_state=seed)
    Xp = pca.fit_transform(Xs)
    dgms_pca = ripser(Xp, maxdim=1)["dgms"]
    pca_stat = max_h1_persistence(dgms_pca)
    return raw_stat, pca_stat


def permute_columns_vectorized(X, rng):
    n, g = X.shape
    rand_vals = rng.random_sample((n, g))
    idx = np.argsort(rand_vals, axis=0)
    return np.take_along_axis(X, idx, axis=0)


def run_gaussian_null(n_samples, n_genes, n_draws, seed=42, n_pca=50):
    gaussian_raw = np.zeros(n_draws)
    gaussian_pca = np.zeros(n_draws)
    t0 = time.time()
    rng_master = np.random.RandomState(seed)
    for i in range(n_draws):
        seed_i = rng_master.randint(0, 2**31 - 1)
        rng_i = np.random.RandomState(seed_i)
        Xg = rng_i.normal(size=(n_samples, n_genes))
        r, p = pipeline_max_h1(Xg, n_pca=n_pca, seed=42)
        gaussian_raw[i] = r
        gaussian_pca[i] = p
        if (i + 1) % max(1, n_draws // 10) == 0:
            print(f"  gaussian null {i+1}/{n_draws} done, elapsed {time.time()-t0:.1f}s")
    print("Gaussian null: raw mean/sd:", gaussian_raw.mean(), gaussian_raw.std())
    print("Gaussian null: pca mean/sd:", gaussian_pca.mean(), gaussian_pca.std())
    return gaussian_raw, gaussian_pca


def run_permutation_null(X_nz, n_perm, seed=42, n_pca=50, n_hvg=2000):
    perm_raw = np.zeros(n_perm)
    perm_pca = np.zeros(n_perm)
    t0 = time.time()
    rng_master = np.random.RandomState(seed)
    for i in range(n_perm):
        seed_i = rng_master.randint(0, 2**31 - 1)
        rng_i = np.random.RandomState(seed_i)
        X_p = permute_columns_vectorized(X_nz, rng_i)
        var_p = X_p.var(axis=0)
        order_p = np.argsort(var_p)
        top_idx_p = order_p[-n_hvg:]
        X_p_hvg = X_p[:, top_idx_p]
        r, p = pipeline_max_h1(X_p_hvg, n_pca=n_pca, seed=42)
        perm_raw[i] = r
        perm_pca[i] = p
        if (i + 1) % max(1, n_perm // 10) == 0:
            print(f"  permutation null {i+1}/{n_perm} done, elapsed {time.time()-t0:.1f}s")
    print("Pipeline-symmetric null: raw mean/sd:", perm_raw.mean(), perm_raw.std())
    print("Pipeline-symmetric null: pca mean/sd:", perm_pca.mean(), perm_pca.std())
    return perm_raw, perm_pca


def estimate_intrinsic_dim_mle(X, k=10):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    dists, _ = nn.kneighbors(X)
    dists = dists[:, 1:]
    log_ratios = np.log(dists[:, -1:] / np.maximum(dists[:, :-1], 1e-15))
    d_hat = (k - 1) / np.sum(log_ratios, axis=1)
    return float(np.median(d_hat))


def compute_p_and_z(observed, null_dist):
    null_mean = null_dist.mean()
    null_sd = null_dist.std()
    z = (observed - null_mean) / null_sd if null_sd > 0 else np.nan
    p = (np.sum(null_dist >= observed) + 1) / (len(null_dist) + 1)
    return null_mean, null_sd, z, p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-gauss", type=int, default=500, help="Gaussian null draws (pre-registered: 500)")
    ap.add_argument("--n-perm", type=int, default=2000, help="Permutation null draws (pre-registered: 2000)")
    ap.add_argument("--skip-assert", action="store_true", help="Print but do not assert against expected values")
    args = ap.parse_args()

    reduced = args.n_gauss != 500 or args.n_perm != 2000
    if reduced:
        print(f"NOTE: running with reduced draw counts (n_gauss={args.n_gauss}, n_perm={args.n_perm}) "
              f"instead of the pre-registered (500, 2000). This is a smoke test, not the confirmatory result.")

    fetch_data()
    X_nz, X_hvg, X_bottom, labels = load_and_preprocess()

    # --- Real-data persistent homology ---
    X_hvg_std = standardize(X_hvg)
    X_bottom_std = standardize(X_bottom)
    pca_hvg = PCA(n_components=50, random_state=42)
    X_hvg_pca = pca_hvg.fit_transform(X_hvg_std)
    print("HVG PCA50 explained var ratio sum:", pca_hvg.explained_variance_ratio_.sum())

    spaces = {"raw_HVG": X_hvg_std, "bottom_variance": X_bottom_std, "PCA50": X_hvg_pca}
    real_diagrams = {name: compute_ph(X, maxdim=2) for name, X in spaces.items()}
    observed_stats = {name: max_h1_persistence(dgms) for name, dgms in real_diagrams.items()}
    print("Observed max-H1 persistence:", observed_stats)

    # --- Intrinsic dimension gate ---
    d_mle_pca50 = estimate_intrinsic_dim_mle(X_hvg_pca, k=10)
    print(f"Intrinsic dimension (MLE, k=10), PCA50 space: {d_mle_pca50:.3f}")

    # --- Null models ---
    n_samples, n_genes = X_hvg.shape
    gaussian_raw, gaussian_pca = run_gaussian_null(n_samples, n_genes, args.n_gauss, seed=42, n_pca=50)
    perm_raw, perm_pca = run_permutation_null(X_nz, args.n_perm, seed=42, n_pca=50, n_hvg=n_genes)

    # --- Statistics table ---
    comparisons = [
        ("Real raw HVG", observed_stats["raw_HVG"], "Gaussian null (raw)", gaussian_raw),
        ("Real raw HVG", observed_stats["raw_HVG"], "Pipeline-symmetric null (raw)", perm_raw),
        ("Real PCA50", observed_stats["PCA50"], "Gaussian null (PCA50)", gaussian_pca),
        ("Real PCA50", observed_stats["PCA50"], "Pipeline-symmetric null (PCA50)", perm_pca),
        ("Real bottom-variance genes", observed_stats["bottom_variance"], "Gaussian null (raw)", gaussian_raw),
        ("Real bottom-variance genes", observed_stats["bottom_variance"], "Pipeline-symmetric null (raw)", perm_raw),
    ]
    rows = []
    for cond_name, obs, null_name, null_dist in comparisons:
        null_mean, null_sd, z, p = compute_p_and_z(obs, null_dist)
        rows.append({
            "condition": cond_name, "observed_max_H1_persistence": obs, "null_model": null_name,
            "null_mean": null_mean, "null_sd": null_sd, "z_score": z, "raw_p_value": p,
            "validation_method": (f"Monte Carlo Gaussian null, n={len(null_dist)}" if "Gaussian" in null_name
                                   else f"Permutation test (per-gene shuffle) + identical pipeline, n={len(null_dist)}"),
        })
    results_table = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(results_table.to_string(index=False))

    # --- Tumor/normal CV-AUC confound check (PCA50 features) ---
    y = (labels == "tumor").astype(int)
    clf = LogisticRegression(max_iter=5000, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = cross_val_score(clf, X_hvg_pca, y, cv=cv, scoring="roc_auc")
    cv_auc_mean = aucs.mean()
    print("5-fold CV AUC (PCA50):", aucs, "mean:", cv_auc_mean, "sd:", aucs.std())

    out_csv = os.path.join(_HERE, "final_results_table_GSE146889.csv")
    results_table.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # --- Verification against the confirmatory-test report ---
    # NOTE: real_max_h1_raw_hvg, real_max_h1_pca50, and d_int_mle_pca50 are computed
    # entirely from the real data (fixed GEO download + fixed seeds) and do NOT depend
    # on n_gauss/n_perm at all -- they are asserted at full tolerance even in a
    # reduced-draw smoke test. Only the null-distribution z-scores/p-values in
    # final_results_table_GSE146889.csv are noisier / not assertable under --n-gauss/--n-perm.
    print("\n=== Verification against report_GSE146889.md ===")
    checks = [
        ("real_max_h1_raw_hvg", observed_stats["raw_HVG"], TOL),
        ("real_max_h1_pca50", observed_stats["PCA50"], TOL),
        ("cv_auc_mean", cv_auc_mean, TOL_CV_AUC),
        ("d_int_mle_pca50", d_mle_pca50, TOL_DIM),
    ]
    all_ok = True
    for name, val, tol in checks:
        exp = EXPECTED[name]
        ok = abs(val - exp) < tol
        all_ok = all_ok and ok
        status = "PASS" if ok else "FAIL"
        print(f"  {name}: got {val:.4f}, expected {exp:.4f}, tol={tol} -> {status}")

    if not args.skip_assert:
        for name, val, tol in checks:
            exp = EXPECTED[name]
            assert abs(val - exp) < tol, f"{name} mismatch: got {val}, expected {exp} (tol {tol})"
        print("\nAll assertions PASSED.")
    if reduced:
        print("\nNOTE: this was a reduced-draw-count run (see NOTE above) -- the real-data numbers "
              "above are exact/assertable regardless, but the null-distribution z-scores/p-values in "
              "the results table are noisier than the pre-registered n_gauss=500/n_perm=2000 run.")


if __name__ == "__main__":
    main()
