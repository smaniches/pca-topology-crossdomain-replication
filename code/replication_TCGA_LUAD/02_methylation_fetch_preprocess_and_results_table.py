#!/usr/bin/env python3
"""
Reproduces the TCGA-LUAD confirmatory replication cohort, methylation layer,
end-to-end: fetch (or reuse a cached copy of) Illumina 450K methylation
beta-values from the GDC public API, preprocess (logit transform),
compute persistent homology (real data, raw top-2000-HVP space AND the
PCA(35) space under BOTH a Euclidean metric [reported for transparency
only] and the Section-5-mandated spectral/graph-Laplacian metric [the
counted primary-family metric, since this cohort's PCA space falls in the
"transitional" intrinsic-dimension regime]), run both null models, the
intrinsic-dimension gate, and write tcga_luad_methylation_results_table.csv.

Cohort: TCGA-LUAD (lung adenocarcinoma), DNA methylation (Illumina Human
Methylation 450K, Beta Value), 36 samples (18 paired tumor + 18 paired
normal -- restricted to the cases that ALSO have paired RNA-seq, matching
the RNA-seq layer's case set). This is part of the 3rd pre-registered
confirmatory test in the TOPOLOGICA PCA-topology program (see
../../prereg/PREREGISTRATION.md).

DISCLOSED DEVIATIONS (see printed notes at runtime and
tcga_luad_report.md Sections 4/5/8 for full detail):
  1. Imputation/logit domain collision: beta-values are bounded in [0,1];
     the blanket "missing -> 0" rule collides with logit(beta) diverging to
     -inf at beta=0. DEVIATION: probes NaN in ALL samples are dropped before
     HVG selection; remaining NaNs are imputed to 0 in beta-space, then all
     values are clipped to [1e-6, 1-1e-6] before the logit transform.
  2. PCA(50) infeasible at n=36 samples (sklearn full-SVD cap =
     min(n_samples, n_features)). DEVIATION: n_components = min(50,
     n_samples-1) = 35 for methylation only.
  3. Section-5 metric gate: this cohort's PCA(35) space falls in the
     "transitional" intrinsic-dimension regime, which the pre-registered
     protocol's Section 5 says requires the spectral (graph-Laplacian
     eigenmap) distance metric rather than raw Euclidean distance for the
     PRIMARY comparison. Euclidean-PCA numbers are still computed and
     reported for transparency, but are NOT part of the primary test family.

Expected results (tolerance 1e-3 unless noted; see
results/replication_TCGA_LUAD/tcga_luad_report.md for the full table):
    real max-H1 persistence, PCA(35, spectral metric) = 0.0 EXACTLY
        -- this is the reported null-space collapse (a real finding, not a
        bug): the spectral embedding of this particular 36-sample point
        cloud produces zero H1 persistence.

Usage:
    python3 02_methylation_fetch_preprocess_and_results_table.py
        # full run: N_GAUSS=500, N_PERM=2000 (as pre-registered).
        # methylation null draws include an eigendecomposition per draw
        # (spectral distance) on top of the ripser call -- slower per draw
        # than the RNA-seq layer; ~30-60 min on a single core for the full run.

    python3 02_methylation_fetch_preprocess_and_results_table.py --n-gauss 15 --n-perm 15
        # fast smoke test (a few minutes) -- NOT the pre-registered draw
        # count; use only to verify the pipeline runs end-to-end.

Data source: by default this script looks for a cached copy of the beta
matrix + sample metadata under ./data/ (see --data-dir: expects
meth_beta_matrix.pkl [probes x samples, pandas DataFrame pickle] and
meth_sample_meta.csv). If absent, it re-fetches directly from the public GDC
API (documented explicitly, not assumed) -- same manifest/download mechanics
as 01_rnaseq_fetch_preprocess_and_results_table.py but filtered to
data_category="DNA Methylation", data_type="Methylation Beta Value",
platform="Illumina Human Methylation 450", restricted to the case overlap
with the RNA-seq layer's paired tumor+normal cases.
"""
import argparse
import os
import sys
import time
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (select_hvg, standardize_samples_x_features, dim_gate,
                      max_h1_persistence, max_h1_persistence_from_D,
                      compute_spectral_distance, compute_z_and_p, benjamini_hochberg)

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HERE = os.path.dirname(os.path.abspath(__file__))
GDC_API = "https://api.gdc.cancer.gov"

EXPECTED = {"real_max_h1_pca35_spectral": 0.0}
TOL = 1e-9  # exact null-space collapse; near-zero tolerance


def fetch_from_gdc(data_dir, rna_case_overlap=None):
    """Re-fetch Illumina 450K methylation beta-values for TCGA-LUAD cases
    with both tumor and matched-normal samples, directly from the public
    GDC API, restricted to the case overlap with the RNA-seq layer if
    provided (matching the original replication's case selection)."""
    import requests

    print("Querying GDC for TCGA-LUAD methylation (Illumina 450K) file manifest...")
    filt = {
        "op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["DNA Methylation"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
            {"op": "in", "content": {"field": "platform", "value": ["Illumina Human Methylation 450"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                      "value": ["Primary Tumor", "Solid Tissue Normal"]}},
        ]
    }
    r = requests.post(f"{GDC_API}/files", json={
        "filters": filt,
        "fields": "file_id,file_name,file_size,cases.submitter_id,cases.samples.sample_type",
        "size": 800, "format": "json",
    }, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print("Total methylation files:", len(hits))

    meth_case_files = defaultdict(dict)
    for h in hits:
        case = h["cases"][0]
        sub = case["submitter_id"]
        stype = case["samples"][0]["sample_type"]
        meth_case_files[sub].setdefault(stype, []).append((h["file_id"], h["file_size"]))
    meth_both = {k: v for k, v in meth_case_files.items()
                 if "Primary Tumor" in v and "Solid Tissue Normal" in v}
    print("Methylation cases with both tumor+normal:", len(meth_both))

    case_set = set(meth_both.keys())
    if rna_case_overlap is not None:
        case_set = case_set & set(rna_case_overlap)
        print("Restricted to overlap with RNA-seq case set:", len(case_set))

    meth_files = []
    for case in sorted(case_set):
        d = meth_both[case]
        for stype in ["Primary Tumor", "Solid Tissue Normal"]:
            fid, fsize = d[stype][0]
            meth_files.append((fid, case, stype))
    print("Methylation files to fetch:", len(meth_files))

    meth_dir = os.path.join(data_dir, "_raw_meth_files")
    os.makedirs(meth_dir, exist_ok=True)
    for fid, case, stype in meth_files:
        out_path = os.path.join(meth_dir, fid + ".txt")
        if os.path.isfile(out_path):
            continue
        resp = requests.get(f"{GDC_API}/data/{fid}", timeout=120)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            f.write(resp.content)

    meth_cols = {}
    meth_meta = []
    for fid, case, stype in meth_files:
        path = os.path.join(meth_dir, fid + ".txt")
        s = pd.read_csv(path, sep="\t", header=None, index_col=0, names=["probe", "beta"])["beta"]
        col_name = f"{case}__{stype.replace(' ', '_')}"
        meth_cols[col_name] = s
        meth_meta.append({"sample": col_name, "case": case, "sample_type": stype, "file_id": fid})

    meth_df = pd.DataFrame(meth_cols)
    meth_meta_df = pd.DataFrame(meth_meta)
    print("Methylation matrix (probes x samples):", meth_df.shape)
    print(meth_meta_df["sample_type"].value_counts())
    return meth_df, meth_meta_df


def load_data(data_dir):
    pkl_path = os.path.join(data_dir, "meth_beta_matrix.pkl")
    meta_path = os.path.join(data_dir, "meth_sample_meta.csv")
    if os.path.isfile(pkl_path) and os.path.isfile(meta_path):
        print(f"Using cached methylation data from {data_dir}")
        meth_df = pd.read_pickle(pkl_path)
        meta_df = pd.read_csv(meta_path)
        print("Methylation matrix (probes x samples):", meth_df.shape)
        print(meta_df["sample_type"].value_counts())
        return meth_df, meta_df
    print(f"No cached data found under {data_dir}; fetching fresh from the GDC public API.")
    os.makedirs(data_dir, exist_ok=True)
    # Try to restrict to the RNA-seq layer's case set if that data is cached alongside.
    rna_meta_path = os.path.join(data_dir, "rna_sample_meta.csv")
    rna_case_overlap = None
    if os.path.isfile(rna_meta_path):
        rna_case_overlap = pd.read_csv(rna_meta_path)["case"].unique().tolist()
    meth_df, meta_df = fetch_from_gdc(data_dir, rna_case_overlap=rna_case_overlap)
    meth_df.to_pickle(pkl_path)
    meta_df.to_csv(meta_path, index=False)
    return meth_df, meta_df


def run_gaussian_null(n_samples, n_features, n_draws, n_pca, seed=42):
    raw_vals, pca_euc_vals, pca_spec_vals = [], [], []
    t0 = time.time()
    rng = np.random.RandomState(seed)
    for i in range(n_draws):
        Xg = rng.standard_normal(size=(n_samples, n_features))
        mu = Xg.mean(axis=0, keepdims=True)
        sd = Xg.std(axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        Xg = (Xg - mu) / sd
        raw_vals.append(max_h1_persistence(Xg)["max_H1"])
        Xg_pca = PCA(n_components=n_pca, random_state=42).fit_transform(Xg)
        pca_euc_vals.append(max_h1_persistence(Xg_pca)["max_H1"])
        D_spec = compute_spectral_distance(Xg_pca)
        pca_spec_vals.append(max_h1_persistence_from_D(D_spec)["max_H1"])
        if (i + 1) % max(1, n_draws // 10) == 0:
            print(f"  gaussian null {i+1}/{n_draws} done, elapsed {time.time()-t0:.1f}s")
    return np.array(raw_vals), np.array(pca_euc_vals), np.array(pca_spec_vals)


def run_permutation_null(mat_probes_x_samples, n_perm, n_pca, seed=42):
    raw_vals, pca_euc_vals, pca_spec_vals = [], [], []
    t0 = time.time()
    rng = np.random.RandomState(seed)
    arr = mat_probes_x_samples.values.copy()
    for i in range(n_perm):
        perm_arr = arr.copy()
        for row in range(perm_arr.shape[0]):
            rng.shuffle(perm_arr[row])
        Xp = perm_arr.T
        mu = Xp.mean(axis=0, keepdims=True)
        sd = Xp.std(axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        Xp = (Xp - mu) / sd
        raw_vals.append(max_h1_persistence(Xp)["max_H1"])
        Xp_pca = PCA(n_components=n_pca, random_state=42).fit_transform(Xp)
        pca_euc_vals.append(max_h1_persistence(Xp_pca)["max_H1"])
        D_spec = compute_spectral_distance(Xp_pca)
        pca_spec_vals.append(max_h1_persistence_from_D(D_spec)["max_H1"])
        if (i + 1) % max(1, n_perm // 10) == 0:
            print(f"  permutation null {i+1}/{n_perm} done, elapsed {time.time()-t0:.1f}s")
    return np.array(raw_vals), np.array(pca_euc_vals), np.array(pca_spec_vals)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-gauss", type=int, default=500)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--data-dir", type=str, default=os.path.join(_HERE, "data"))
    ap.add_argument("--skip-assert", action="store_true")
    args = ap.parse_args()

    reduced = args.n_gauss != 500 or args.n_perm != 2000
    if reduced:
        print(f"NOTE: running with reduced draw counts (n_gauss={args.n_gauss}, n_perm={args.n_perm}) "
              f"instead of the pre-registered (500, 2000). This is a smoke test, not the confirmatory result.")

    meth_df, meta_df = load_data(args.data_dir)

    print("\n--- DISCLOSED DEVIATION 1: imputation/logit domain collision (see docstring) ---")
    nan_frac = meth_df.isna().mean().mean()
    n_allnan_rows = meth_df.isna().all(axis=1).sum()
    print(f"NaN fraction: {nan_frac:.4f} | Probes NaN in ALL samples: {n_allnan_rows} / {len(meth_df)}")
    meth_df_clean = meth_df.loc[~meth_df.isna().all(axis=1)]
    print("After dropping all-NaN probes:", meth_df_clean.shape)
    meth_imputed = meth_df_clean.fillna(0.0)
    eps = 1e-6
    meth_clipped = meth_imputed.clip(lower=eps, upper=1 - eps)
    logit_meth = np.log(meth_clipped / (1 - meth_clipped))
    print("logit_meth shape:", logit_meth.shape,
          "| any inf/nan:", np.isinf(logit_meth.values).sum(), np.isnan(logit_meth.values).sum())

    var_meth = logit_meth.var(axis=1)
    n_zero_var_meth = (var_meth == 0).sum()
    print("Zero-variance probes:", n_zero_var_meth, "/", logit_meth.shape[0])

    top_hvp = select_hvg(logit_meth, n_top=2000, bottom=False)
    nonzero_var_probes = logit_meth.loc[var_meth > 0]
    bottom_vp_dev = select_hvg(nonzero_var_probes, n_top=2000, bottom=True)

    X_meth_raw = standardize_samples_x_features(top_hvp)
    X_meth_bottom = standardize_samples_x_features(bottom_vp_dev)

    print("\n--- DISCLOSED DEVIATION 2: PCA(35) not PCA(50), n_samples-1 cap (see docstring) ---")
    meth_pca_k = min(50, X_meth_raw.shape[0] - 1)
    print("Methylation PCA components used:", meth_pca_k)
    pca_meth = PCA(n_components=meth_pca_k, random_state=42)
    X_meth_pca = pca_meth.fit_transform(X_meth_raw)
    print("X_meth_pca:", X_meth_pca.shape, "explained var:", pca_meth.explained_variance_ratio_.sum())

    print("\n--- Intrinsic dimension gate (Methylation) ---")
    for name, X in [("raw_hvp", X_meth_raw), ("bottom_var", X_meth_bottom), (f"pca{meth_pca_k}", X_meth_pca)]:
        print(name, dim_gate(X))

    print("\n--- DISCLOSED DEVIATION 3: Section-5 spectral-metric gate on PCA space (see docstring) ---")
    print("Computing real-data persistent homology (raw HVP + PCA35 Euclidean + PCA35 spectral)...")
    real_raw = max_h1_persistence(X_meth_raw)
    real_bottom = max_h1_persistence(X_meth_bottom)
    real_pca_euc = max_h1_persistence(X_meth_pca)
    D_pca_spectral = compute_spectral_distance(X_meth_pca)
    real_pca_spec = max_h1_persistence_from_D(D_pca_spectral)
    print("Real max-H1: raw_hvp =", real_raw["max_H1"], "| bottom_var =", real_bottom["max_H1"],
          "| pca35(euclidean, transparency only) =", real_pca_euc["max_H1"],
          "| pca35(SPECTRAL, Section-5-mandated) =", real_pca_spec["max_H1"])

    print(f"\nRunning Gaussian null (n={args.n_gauss})...")
    gauss_raw, gauss_pca_euc, gauss_pca_spec = run_gaussian_null(
        X_meth_raw.shape[0], 2000, args.n_gauss, n_pca=meth_pca_k, seed=42)

    print(f"\nRunning permutation null, top-HVP probes (n={args.n_perm})...")
    perm_top_raw, perm_top_pca_euc, perm_top_pca_spec = run_permutation_null(
        top_hvp, args.n_perm, n_pca=meth_pca_k, seed=43)

    print(f"\nRunning permutation null, bottom-variance probes (n={args.n_perm})...")
    perm_bot_raw, perm_bot_pca_euc, perm_bot_pca_spec = run_permutation_null(
        bottom_vp_dev, args.n_perm, n_pca=meth_pca_k, seed=44)

    meth_results = {}
    meth_results["raw_hvp__gaussian"] = compute_z_and_p(real_raw["max_H1"], gauss_raw)
    meth_results["raw_hvp__perm"] = compute_z_and_p(real_raw["max_H1"], perm_top_raw)
    meth_results["pca35_euclidean__gaussian"] = compute_z_and_p(real_pca_euc["max_H1"], gauss_pca_euc)
    meth_results["pca35_euclidean__perm"] = compute_z_and_p(real_pca_euc["max_H1"], perm_top_pca_euc)
    meth_results["pca35_spectral__gaussian"] = compute_z_and_p(real_pca_spec["max_H1"], gauss_pca_spec)
    meth_results["pca35_spectral__perm"] = compute_z_and_p(real_pca_spec["max_H1"], perm_top_pca_spec)
    meth_results["bottom_var__gaussian"] = compute_z_and_p(real_bottom["max_H1"], gauss_raw)
    meth_results["bottom_var__perm"] = compute_z_and_p(real_bottom["max_H1"], perm_bot_raw)
    for k, v in meth_results.items():
        print(k, v)

    # Primary family (per Section 5 gate): raw_hvp + pca35_SPECTRAL only (euclidean-PCA
    # is transparency-only, not counted -- see docstring deviation 3).
    primary_names = ["raw_hvp__gaussian", "raw_hvp__perm", "pca35_spectral__gaussian", "pca35_spectral__perm"]
    pvals = [meth_results[n]["p_empirical"] for n in primary_names]
    zscores = [meth_results[n]["z"] for n in primary_names]
    rejected, p_adj = benjamini_hochberg(pvals, alpha=0.05)
    rows = [{"test": n, "real_max_H1": meth_results[n]["real"], "null_mean": meth_results[n]["null_mean"],
             "null_std": meth_results[n]["null_std"], "z_score": z, "p_empirical": p,
             "p_BH": padj, "reject_at_0.05_BH": rej,
             "passes_section7_threshold_z_gt_3": z is not None and z > 3.0}
            for n, z, p, padj, rej in zip(primary_names, zscores, pvals, p_adj, rejected)]
    summary_df = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print("\n--- Primary test family (Section-5-gated: raw + PCA35-SPECTRAL only) ---")
    print(summary_df.to_string(index=False))

    out_csv = os.path.join(_HERE, "tcga_luad_methylation_results_table.csv")
    summary_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # NOTE: real_max_h1_pca35_spectral is computed entirely from the real data
    # (fixed seeds, deterministic PCA + spectral embedding) and does not depend on
    # n_gauss/n_perm.
    print("\n=== Verification against tcga_luad_report.md (Methylation layer) ===")
    checks = [("real_max_h1_pca35_spectral", real_pca_spec["max_H1"], TOL)]
    for name, val, tol in checks:
        exp = EXPECTED[name]
        ok = abs(val - exp) < tol
        print(f"  {name}: got {val:.6f}, expected {exp:.6f}, tol={tol} -> {'PASS' if ok else 'FAIL'} "
              f"(this is the reported null-space collapse under PCA(35)+spectral metric -- a real "
              f"finding, not expected to be nonzero)")
    if not args.skip_assert:
        for name, val, tol in checks:
            exp = EXPECTED[name]
            assert abs(val - exp) < tol, f"{name} mismatch: got {val}, expected {exp} (tol {tol})"
        print("\nAll assertions PASSED.")
    if reduced:
        print("\nNOTE: reduced-draw-count run -- the real-data max-H1 numbers above are "
              "exact/assertable regardless, but null-distribution z-scores/p-values are noisier "
              "than the pre-registered n_gauss=500/n_perm=2000 run.")


if __name__ == "__main__":
    main()
