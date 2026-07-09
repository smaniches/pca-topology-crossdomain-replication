#!/usr/bin/env python3
"""
Reproduces the CPTAC-CCRCC confirmatory replication cohort end-to-end:
fetch (or reuse a cached copy of) the public PDC log2-ratio proteomics
matrix, preprocess, compute persistent homology (real data, raw top-2000-HVG
and PCA(50) spaces) via gudhi's Vietoris-Rips implementation, run both null
models (Gaussian and pipeline-symmetric permutation), the intrinsic-dimension
gate, and the tumor/normal CV-AUC confound check, and write
cptac_final_results_table.csv.

Cohort: CPTAC-CCRCC (clear cell renal cell carcinoma), proteomics,
PDC study PDC000127, 194 samples (110 Primary Tumor, 84 Solid Tissue Normal).
This is the 1st pre-registered confirmatory test in the TOPOLOGICA
PCA-topology program (see ../../prereg/PREREGISTRATION.md for the locked
protocol).

DISCLOSED DEVIATION (documented explicitly, per the pre-registration's own
instruction for real-data-quirk deviations -- see Section 4/8 of the
pre-registration and the printed note below): the PDC public API only
exposes CCRCC protein abundance as log2(ratio-to-common-reference), not raw
linear intensities, so the pre-registered "median-normalize then log2" rule
is applied as median-CENTERING in log-space (the log2 step is skipped since
the data already arrives log-transformed).

Expected results (tolerance 1e-3 unless noted; see
results/replication_CPTAC_CCRCC/cptac_ccrcc_report.md for the full table):
    real max-H1 persistence, raw top-2000-HVG   = 3.833
    real max-H1 persistence, PCA(50)            = 4.683
    5-fold CV AUC (tumor vs normal, PCA50)       = 1.000

Usage:
    python3 01_fetch_preprocess_and_results_table.py
        # full run: N_GAUSS=500, N_PERM=2000 (as pre-registered).
        # gudhi's Rips-complex persistent homology on ~194-sample point
        # clouds is fast per draw, but 2500 total draws still takes on the
        # order of 20-40 minutes on a single core.

    python3 01_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
        # fast smoke test (a few minutes) -- NOT the pre-registered draw
        # count; use only to verify the pipeline runs end-to-end.

Data source: by default this script looks for a cached copy of the
preprocessed log2-ratio matrix under ./data/ (see --data-dir). If absent, it
re-fetches directly from the public Proteomic Data Commons (PDC) GraphQL API
(documented explicitly, not assumed):
    https://pdc.cancer.gov/graphql
    query: quantDataMatrix(pdc_study_id: "PDC000127", data_type: "log2_ratio", acceptDUA: true)
and the matching biospecimen metadata (biospecimenPerStudy) to label samples
Primary Tumor / Solid Tissue Normal.
"""
import argparse
import os
import time
import warnings

import numpy as np
import pandas as pd
import gudhi
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HERE = os.path.dirname(os.path.abspath(__file__))

PDC_URL = "https://pdc.cancer.gov/graphql"
PDC_STUDY_ID = "PDC000127"

EXPECTED = {
    "real_max_h1_raw": 3.833,
    "real_max_h1_pca50": 4.683,
    "cv_auc_mean": 1.000,
}
TOL = 1e-3
TOL_CV_AUC = 0.02


def fetch_from_pdc():
    """Re-fetch the log2-ratio quant matrix and biospecimen metadata directly
    from the public PDC GraphQL API. Used only if no cached data/ copy is
    present (see main())."""
    import requests

    print(f"Fetching quantDataMatrix for {PDC_STUDY_ID} from {PDC_URL} ...")
    quant_query = ('{ quantDataMatrix( pdc_study_id: "' + PDC_STUDY_ID +
                   '" data_type: "log2_ratio" acceptDUA: true ) }')
    rq = requests.post(PDC_URL, json={"query": quant_query}, timeout=180)
    rq.raise_for_status()
    matrix = rq.json()["data"]["quantDataMatrix"]
    df = pd.DataFrame(matrix[1:], columns=matrix[0]).set_index("Gene/Aliquot")
    print("Raw quant matrix (proteins x aliquots):", df.shape)

    q_bio = ('{ biospecimenPerStudy(pdc_study_id: "' + PDC_STUDY_ID + '", acceptDUA: true) '
             '{ aliquot_id sample_id case_id sample_submitter_id aliquot_submitter_id sample_type } }')
    rb = requests.post(PDC_URL, json={"query": q_bio}, timeout=120)
    rb.raise_for_status()
    bio = rb.json()["data"]["biospecimenPerStudy"]
    subid_to_type = {b["aliquot_submitter_id"]: b["sample_type"] for b in bio}

    col_subids = [c.split(":")[1] for c in df.columns]
    sample_types = [subid_to_type.get(s, "Unknown") for s in col_subids]

    keep_mask = [st in ("Primary Tumor", "Solid Tissue Normal") for st in sample_types]
    df_use = df.loc[:, keep_mask].copy()
    labels = pd.Series([st for st, k in zip(sample_types, keep_mask) if k], index=df_use.columns)
    df_use = df_use.apply(pd.to_numeric, errors="coerce")

    X_raw = df_use.T  # samples x proteins
    X_raw.index = [c.split(":")[1] for c in X_raw.index]
    labels.index = X_raw.index
    print("Filtered matrix (samples x proteins):", X_raw.shape, "| label counts:", labels.value_counts().to_dict())
    return X_raw, labels


def load_data(data_dir):
    """Load the log2-ratio matrix + labels from a local data/ directory if
    present (recommended: avoids re-hitting the PDC API on every run and
    matches what was used to produce the confirmatory-test report), else
    fetch fresh from PDC."""
    raw_path = os.path.join(data_dir, "cptac_ccrcc_log2ratio_raw.csv")
    labels_path = os.path.join(data_dir, "cptac_ccrcc_labels.csv")
    if os.path.isfile(raw_path) and os.path.isfile(labels_path):
        print(f"Using cached preprocessed data from {data_dir}")
        X_raw = pd.read_csv(raw_path, index_col=0)
        labels = pd.read_csv(labels_path, index_col=0).iloc[:, 0]
        labels.index = X_raw.index
        print("Loaded matrix (samples x proteins):", X_raw.shape, "| label counts:", labels.value_counts().to_dict())
        return X_raw, labels
    print(f"No cached data found under {data_dir}; fetching fresh from the PDC public API.")
    os.makedirs(data_dir, exist_ok=True)
    X_raw, labels = fetch_from_pdc()
    X_raw.to_csv(raw_path)
    labels.to_csv(labels_path)
    return X_raw, labels


def preprocess(X_raw):
    print("\n--- DISCLOSED DEVIATION (pre-registration Section 4/8) ---")
    print("PDC exposes only log2(ratio-to-reference), not raw linear intensity; "
          "'median-normalize then log2' is applied as median-CENTERING in log-space, "
          "skipping the redundant log2 step.")
    X_imputed = X_raw.fillna(0.0)  # missing/sentinel -> 0, per the imputation rule
    sample_medians = X_imputed.median(axis=1)
    print("Sample median log-ratio range (pre-centering):", sample_medians.min(), sample_medians.max())
    X_norm = X_imputed.sub(sample_medians, axis=0)
    print("Post-centering sample medians (should be ~0):", X_norm.median(axis=1).abs().max())

    protein_var = X_norm.var(axis=0)
    top2000 = protein_var.sort_values(ascending=False).index[:2000]
    bottom2000 = protein_var.sort_values(ascending=True).index[:2000]
    X_top = X_norm[top2000].values
    X_bottom = X_norm[bottom2000].values
    print("Top-2000-HVG shape:", X_top.shape, "| bottom-2000 shape:", X_bottom.shape)
    return X_norm, X_top, X_bottom


def estimate_intrinsic_dim_mle(X, k=10):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    dists, _ = nn.kneighbors(X)
    dists = dists[:, 1:]
    log_ratios = np.log(dists[:, -1:] / np.maximum(dists[:, :-1], 1e-15))
    d_hat = (k - 1) / np.sum(log_ratios, axis=1)
    return float(np.median(d_hat))


def max_h1_pers_gudhi(X):
    D = squareform(pdist(X, metric="euclidean"))
    mx = D.max()
    if mx <= 0:
        return 0.0
    rips = gudhi.RipsComplex(distance_matrix=D, max_edge_length=mx * 1.001)
    st = rips.create_simplex_tree(max_dimension=2)
    st.compute_persistence()
    h1 = st.persistence_intervals_in_dimension(1)
    h1 = np.array([(b, d) for b, d in h1 if np.isfinite(d)])
    if len(h1) == 0:
        return 0.0
    return float((h1[:, 1] - h1[:, 0]).max())


def run_gaussian_null(n_samples, n_feat, n_draws, seed=42):
    raw_vals = np.zeros(n_draws)
    pca_vals = np.zeros(n_draws)
    t0 = time.time()
    rng = np.random.RandomState(seed)
    for i in range(n_draws):
        G = rng.normal(size=(n_samples, n_feat))
        Gs = StandardScaler().fit_transform(G)
        raw_vals[i] = max_h1_pers_gudhi(Gs)
        Gp = PCA(n_components=50, random_state=42).fit_transform(Gs)
        pca_vals[i] = max_h1_pers_gudhi(Gp)
        if (i + 1) % max(1, n_draws // 10) == 0:
            print(f"  gaussian null {i+1}/{n_draws} done, elapsed {time.time()-t0:.1f}s")
    return raw_vals, pca_vals


def run_permutation_null(X_top_raw, n_perm, seed=42):
    raw_vals = np.zeros(n_perm)
    pca_vals = np.zeros(n_perm)
    t0 = time.time()
    rng = np.random.RandomState(seed)
    for i in range(n_perm):
        Xp = X_top_raw.copy()
        for j in range(Xp.shape[1]):
            Xp[:, j] = rng.permutation(Xp[:, j])
        Xs = StandardScaler().fit_transform(Xp)
        raw_vals[i] = max_h1_pers_gudhi(Xs)
        Xpp = PCA(n_components=50, random_state=42).fit_transform(Xs)
        pca_vals[i] = max_h1_pers_gudhi(Xpp)
        if (i + 1) % max(1, n_perm // 10) == 0:
            print(f"  permutation null {i+1}/{n_perm} done, elapsed {time.time()-t0:.1f}s")
    return raw_vals, pca_vals


def exact_p(observed, null_dist):
    n = len(null_dist)
    count = np.sum(null_dist >= observed)
    return (count + 1) / (n + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-gauss", type=int, default=500, help="Gaussian null draws (pre-registered: 500)")
    ap.add_argument("--n-perm", type=int, default=2000, help="Permutation null draws (pre-registered: 2000)")
    ap.add_argument("--data-dir", type=str, default=os.path.join(_HERE, "data"),
                     help="Directory holding cached cptac_ccrcc_log2ratio_raw.csv / cptac_ccrcc_labels.csv")
    ap.add_argument("--skip-assert", action="store_true")
    args = ap.parse_args()

    reduced = args.n_gauss != 500 or args.n_perm != 2000
    if reduced:
        print(f"NOTE: running with reduced draw counts (n_gauss={args.n_gauss}, n_perm={args.n_perm}) "
              f"instead of the pre-registered (500, 2000). This is a smoke test, not the confirmatory result.")

    X_raw, labels = load_data(args.data_dir)
    X_norm, X_top, X_bottom = preprocess(X_raw)

    X_top_std = StandardScaler().fit_transform(X_top)
    X_bottom_std = StandardScaler().fit_transform(X_bottom)
    pca = PCA(n_components=50, random_state=42)
    X_pca50 = pca.fit_transform(X_top_std)
    print("PCA50 explained var ratio sum:", pca.explained_variance_ratio_.sum())

    d_int_mle_pca = estimate_intrinsic_dim_mle(X_pca50, k=10)
    print(f"Intrinsic dimension (MLE, k=10), PCA50 space: {d_int_mle_pca:.3f}")

    print("\nComputing real-data persistent homology (gudhi Rips complex)...")
    real_raw_stat = max_h1_pers_gudhi(X_top_std)
    real_pca_stat = max_h1_pers_gudhi(X_pca50)
    print("Real max-H1-persistence: raw top-HVG =", real_raw_stat, " PCA50 =", real_pca_stat)

    print(f"\nRunning Gaussian null (n={args.n_gauss})...")
    gauss_raw, gauss_pca = run_gaussian_null(X_top.shape[0], X_top.shape[1], args.n_gauss, seed=42)
    print("Gaussian null: raw mean/sd:", gauss_raw.mean(), gauss_raw.std(),
          "| pca mean/sd:", gauss_pca.mean(), gauss_pca.std())

    print(f"\nRunning pipeline-symmetric permutation null (n={args.n_perm})...")
    perm_raw, perm_pca = run_permutation_null(X_top, args.n_perm, seed=42)
    print("Permutation null: raw mean/sd:", perm_raw.mean(), perm_raw.std(),
          "| pca mean/sd:", perm_pca.mean(), perm_pca.std())

    results = {
        "Real raw-HVG vs Gaussian null": (real_raw_stat, gauss_raw),
        "Real PCA50 vs Gaussian null": (real_pca_stat, gauss_pca),
        "Real raw-HVG vs pipeline-symmetric null": (real_raw_stat, perm_raw),
        "Real PCA50 vs pipeline-symmetric null": (real_pca_stat, perm_pca),
    }
    rows = []
    for label, (obs, null_dist) in results.items():
        p = exact_p(obs, null_dist)
        z = (obs - null_dist.mean()) / null_dist.std() if null_dist.std() > 0 else np.nan
        rows.append({"comparison": label, "observed": obs, "null_mean": null_dist.mean(),
                     "null_sd": null_dist.std(), "z_score": z, "p_value": p, "n_null": len(null_dist)})
        print(f"{label}: observed={obs:.4f}, null_mean={null_dist.mean():.4f}, "
              f"null_sd={null_dist.std():.4f}, z={z:.3f}, p={p:.5f}")

    final_df = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(final_df.to_string(index=False))

    # --- Tumor/normal CV-AUC confound check (PCA50 features) ---
    y = (labels.loc[X_norm.index] == "Primary Tumor").astype(int).values
    clf = LogisticRegression(max_iter=2000, random_state=42)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = cross_val_score(clf, X_pca50, y, cv=skf, scoring="roc_auc")
    cv_auc_mean = aucs.mean()
    print("5-fold CV AUC (PCA50, tumor vs normal):", aucs, "mean:", cv_auc_mean)

    out_csv = os.path.join(_HERE, "cptac_final_results_table.csv")
    final_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # NOTE: real_max_h1_raw and real_max_h1_pca50 are computed entirely from the real
    # data (fixed seeds, PCA is deterministic) and do NOT depend on n_gauss/n_perm --
    # asserted at full tolerance even under a reduced-draw smoke test. Only the
    # null-distribution z-scores/p-values are noisier under reduced draw counts.
    print("\n=== Verification against cptac_ccrcc_report.md ===")
    checks = [
        ("real_max_h1_raw", real_raw_stat, TOL),
        ("real_max_h1_pca50", real_pca_stat, TOL),
        ("cv_auc_mean", cv_auc_mean, TOL_CV_AUC),
    ]
    for name, val, tol in checks:
        exp = EXPECTED[name]
        ok = abs(val - exp) < tol
        status = "PASS" if ok else "FAIL"
        print(f"  {name}: got {val:.4f}, expected {exp:.4f}, tol={tol} -> {status}")

    if not args.skip_assert:
        for name, val, tol in checks:
            exp = EXPECTED[name]
            assert abs(val - exp) < tol, f"{name} mismatch: got {val}, expected {exp} (tol {tol})"
        print("\nAll assertions PASSED.")
    if reduced:
        print("\nNOTE: this was a reduced-draw-count run -- the real-data numbers above are "
              "exact/assertable regardless, but the null-distribution z-scores/p-values in the "
              "results table are noisier than the pre-registered n_gauss=500/n_perm=2000 run.")


if __name__ == "__main__":
    main()
