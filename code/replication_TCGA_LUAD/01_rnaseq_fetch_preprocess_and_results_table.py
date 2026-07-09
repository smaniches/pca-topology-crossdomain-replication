#!/usr/bin/env python3
"""
Reproduces the TCGA-LUAD confirmatory replication cohort, RNA-seq layer,
end-to-end: fetch (or reuse a cached copy of) STAR-Counts FPKM gene
expression from the GDC public API, preprocess, compute persistent
homology (real data, raw top-2000-HVG and PCA(50) spaces), run both null
models, the intrinsic-dimension gate, and the tumor/normal CV-AUC confound
check, and write tcga_luad_rnaseq_results_table.csv.

Cohort: TCGA-LUAD (lung adenocarcinoma), RNA-seq, 116 samples (58 paired
tumor + 58 paired normal, STAR-Counts FPKM-unstranded, Primary Tumor /
Solid Tissue Normal only). This is the 3rd pre-registered confirmatory test
in the TOPOLOGICA PCA-topology program (see ../../prereg/PREREGISTRATION.md).

DISCLOSED DEVIATION (RNA-seq, bottom-variance selection): 5556/60660 genes
have EXACTLY zero variance across all 116 samples (never detected in any
sample), exceeding the pre-registered 2000-gene bottom-variance selection
size. Literal "bottom 2000 by variance" is therefore a degenerate all-zero
point cloud with trivial topology by construction. DEVIATION APPLIED:
bottom-variance genes are selected as the bottom 2000 among genes with
variance > 0 (i.e. detected in at least one sample). This is reported here
explicitly, not silently absorbed.

Expected results (tolerance 1e-3 unless noted; see
results/replication_TCGA_LUAD/tcga_luad_report.md for the full table):
    real max-H1 persistence, raw top-2000-HVG  = 8.327
    real max-H1 persistence, PCA(50)           = 7.433
    5-fold CV AUC (tumor vs normal, PCA50)      = 0.998

Usage:
    python3 01_rnaseq_fetch_preprocess_and_results_table.py
        # full run: N_GAUSS=500, N_PERM=2000 (as pre-registered).
        # ~2500 total ripser calls on ~116x2000 point clouds: on the order
        # of 20-40 minutes on a single core.

    python3 01_rnaseq_fetch_preprocess_and_results_table.py --n-gauss 20 --n-perm 20
        # fast smoke test (a few minutes) -- NOT the pre-registered draw
        # count; use only to verify the pipeline runs end-to-end.

Data source: by default this script looks for a cached copy of the raw
FPKM matrix + sample metadata under ./data/ (see --data-dir: expects
rna_fpkm_matrix.pkl [genes x samples, pandas DataFrame pickle] and
rna_sample_meta.csv). If absent, it re-fetches directly from the public GDC
API (documented explicitly, not assumed):
    https://api.gdc.cancer.gov/files   (query: TCGA-LUAD, STAR - Counts,
        Gene Expression Quantification, Primary Tumor / Solid Tissue Normal,
        restricted to the 58 cases with BOTH tumor and matched-normal samples)
    https://api.gdc.cancer.gov/data/{file_id}   (per-file download)
"""
import argparse
import json
import os
import subprocess
import sys
import time
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (select_hvg, standardize_samples_x_features, dim_gate,
                      max_h1_persistence, compute_z_and_p, benjamini_hochberg)

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HERE = os.path.dirname(os.path.abspath(__file__))
GDC_API = "https://api.gdc.cancer.gov"

EXPECTED = {"real_max_h1_raw_hvg": 8.327, "real_max_h1_pca50": 7.433, "cv_auc_mean": 0.998}
TOL = 1e-3
TOL_CV_AUC = 0.02


def fetch_from_gdc(data_dir):
    """Re-fetch STAR-Counts FPKM RNA-seq for the 58 TCGA-LUAD cases with
    BOTH a Primary Tumor and a matched Solid Tissue Normal sample, directly
    from the public GDC API."""
    import requests

    print("Querying GDC for TCGA-LUAD RNA-seq (STAR - Counts) file manifest...")
    filt = {
        "op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                      "value": ["Primary Tumor", "Solid Tissue Normal"]}},
        ]
    }
    r = requests.post(f"{GDC_API}/files", json={
        "filters": filt,
        "fields": "file_id,file_name,cases.case_id,cases.submitter_id,"
                  "cases.samples.sample_type,cases.samples.submitter_id",
        "size": 700, "format": "json",
    }, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print("Total files:", len(hits))

    case_files = defaultdict(dict)
    for h in hits:
        case = h["cases"][0]
        sub = case["submitter_id"]
        stype = case["samples"][0]["sample_type"]
        case_files[sub].setdefault(stype, []).append(h["file_id"])
    both = {k: v for k, v in case_files.items() if "Primary Tumor" in v and "Solid Tissue Normal" in v}
    print("Cases with both tumor+normal:", len(both))

    rna_files = []
    for case, d in both.items():
        for stype in ["Primary Tumor", "Solid Tissue Normal"]:
            rna_files.append((d[stype][0], case, stype))
    print("RNA-seq files to fetch:", len(rna_files))

    rna_dir = os.path.join(data_dir, "_raw_rna_files")
    os.makedirs(rna_dir, exist_ok=True)
    for fid, case, stype in rna_files:
        out_path = os.path.join(rna_dir, fid + ".tsv")
        if os.path.isfile(out_path):
            continue
        resp = requests.get(f"{GDC_API}/data/{fid}", timeout=120)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            f.write(resp.content)

    samples_meta = []
    expr_cols = {}
    for fid, case, stype in rna_files:
        path = os.path.join(rna_dir, fid + ".tsv")
        df = pd.read_csv(path, sep="\t", skiprows=1)
        df = df[df["gene_id"].str.startswith("ENSG")]
        col_name = f"{case}__{stype.replace(' ', '_')}"
        expr_cols[col_name] = df.set_index("gene_id")["fpkm_unstranded"]
        samples_meta.append({"sample": col_name, "case": case, "sample_type": stype, "file_id": fid})

    expr_df = pd.DataFrame(expr_cols)
    meta_df = pd.DataFrame(samples_meta)
    print("Expression matrix (genes x samples):", expr_df.shape)
    print(meta_df["sample_type"].value_counts())
    return expr_df, meta_df


def load_data(data_dir):
    pkl_path = os.path.join(data_dir, "rna_fpkm_matrix.pkl")
    meta_path = os.path.join(data_dir, "rna_sample_meta.csv")
    if os.path.isfile(pkl_path) and os.path.isfile(meta_path):
        print(f"Using cached RNA-seq data from {data_dir}")
        expr_df = pd.read_pickle(pkl_path)
        meta_df = pd.read_csv(meta_path)
        print("Expression matrix (genes x samples):", expr_df.shape)
        print(meta_df["sample_type"].value_counts())
        return expr_df, meta_df
    print(f"No cached data found under {data_dir}; fetching fresh from the GDC public API.")
    os.makedirs(data_dir, exist_ok=True)
    expr_df, meta_df = fetch_from_gdc(data_dir)
    expr_df.to_pickle(pkl_path)
    meta_df.to_csv(meta_path, index=False)
    return expr_df, meta_df


def run_gaussian_null(n_samples, n_features, n_draws, n_pca, seed=42):
    raw_vals, pca_vals = [], []
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
        pca_vals.append(max_h1_persistence(Xg_pca)["max_H1"])
        if (i + 1) % max(1, n_draws // 10) == 0:
            print(f"  gaussian null {i+1}/{n_draws} done, elapsed {time.time()-t0:.1f}s")
    return np.array(raw_vals), np.array(pca_vals)


def run_permutation_null(mat_genes_x_samples, n_perm, n_pca, seed=42):
    """mat_genes_x_samples: genes x samples DataFrame, already HVG-subset."""
    raw_vals, pca_vals = [], []
    t0 = time.time()
    rng = np.random.RandomState(seed)
    arr = mat_genes_x_samples.values.copy()
    for i in range(n_perm):
        perm_arr = arr.copy()
        for row in range(perm_arr.shape[0]):
            rng.shuffle(perm_arr[row])
        Xp = perm_arr.T  # samples x genes
        mu = Xp.mean(axis=0, keepdims=True)
        sd = Xp.std(axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        Xp = (Xp - mu) / sd
        raw_vals.append(max_h1_persistence(Xp)["max_H1"])
        Xp_pca = PCA(n_components=n_pca, random_state=42).fit_transform(Xp)
        pca_vals.append(max_h1_persistence(Xp_pca)["max_H1"])
        if (i + 1) % max(1, n_perm // 10) == 0:
            print(f"  permutation null {i+1}/{n_perm} done, elapsed {time.time()-t0:.1f}s")
    return np.array(raw_vals), np.array(pca_vals)


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

    expr_df, meta_df = load_data(args.data_dir)

    expr_df = expr_df.fillna(0.0)  # missing/sentinel -> 0
    log_expr = np.log1p(expr_df)   # RNA-seq preprocessing: log1p(FPKM)
    print("log_expr shape (genes x samples):", log_expr.shape)

    top_hvg = select_hvg(log_expr, n_top=2000, bottom=False)

    var_all = log_expr.var(axis=1)
    n_zero_var = (var_all == 0).sum()
    print(f"Genes with exactly zero variance across all {log_expr.shape[1]} samples: "
          f"{n_zero_var} / {log_expr.shape[0]}")
    print("DISCLOSED DEVIATION: bottom-variance genes selected among variance>0 pool only (see docstring).")
    nonzero_var_genes = log_expr.loc[var_all > 0]
    bottom_var_dev = select_hvg(nonzero_var_genes, n_top=2000, bottom=True)

    X_raw_hvg = standardize_samples_x_features(top_hvg)      # 116 x 2000
    X_bottom = standardize_samples_x_features(bottom_var_dev)
    pca = PCA(n_components=50, random_state=42)
    X_pca50 = pca.fit_transform(X_raw_hvg)
    print("X_raw_hvg:", X_raw_hvg.shape, "| X_bottom:", X_bottom.shape,
          "| X_pca50:", X_pca50.shape, "explained var:", pca.explained_variance_ratio_.sum())

    print("\n--- Intrinsic dimension gate (RNA-seq) ---")
    for name, X in [("raw_hvg", X_raw_hvg), ("bottom_var", X_bottom), ("pca50", X_pca50)]:
        print(name, dim_gate(X))

    print("\nComputing real-data persistent homology...")
    real_raw = max_h1_persistence(X_raw_hvg)
    real_bottom = max_h1_persistence(X_bottom)
    real_pca = max_h1_persistence(X_pca50)
    real_raw_h1 = real_raw["max_H1"]
    real_pca_h1 = real_pca["max_H1"]
    print("Real max-H1 persistence: raw_hvg =", real_raw_h1, "| bottom_var =", real_bottom["max_H1"],
          "| pca50 =", real_pca_h1)

    print(f"\nRunning Gaussian null (n={args.n_gauss})...")
    gauss_raw, gauss_pca = run_gaussian_null(X_raw_hvg.shape[0], 2000, args.n_gauss, n_pca=50, seed=42)

    print(f"\nRunning permutation null, top-HVG genes (n={args.n_perm})...")
    perm_top_raw, perm_top_pca = run_permutation_null(top_hvg, args.n_perm, n_pca=50, seed=43)

    print(f"\nRunning permutation null, bottom-variance genes (n={args.n_perm})...")
    perm_bot_raw, perm_bot_pca = run_permutation_null(bottom_var_dev, args.n_perm, n_pca=50, seed=44)

    rna_results = {}
    rna_results["raw_hvg__gaussian"] = compute_z_and_p(real_raw_h1, gauss_raw)
    rna_results["raw_hvg__perm"] = compute_z_and_p(real_raw_h1, perm_top_raw)
    rna_results["pca50__gaussian"] = compute_z_and_p(real_pca_h1, gauss_pca)
    rna_results["pca50__perm"] = compute_z_and_p(real_pca_h1, perm_top_pca)
    rna_results["bottom_var__gaussian"] = compute_z_and_p(real_bottom["max_H1"], gauss_raw)
    rna_results["bottom_var__perm"] = compute_z_and_p(real_bottom["max_H1"], perm_bot_raw)
    for k, v in rna_results.items():
        print(k, v)

    # --- Tumor/normal CV-AUC confound check (PCA50 features) ---
    rna_labels = meta_df["sample_type"].values
    y = np.array([1 if l == "Primary Tumor" else 0 for l in rna_labels])
    clf = LogisticRegression(max_iter=2000, random_state=42)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = cross_val_score(clf, X_pca50, y, cv=skf, scoring="roc_auc")
    cv_auc_mean = aucs.mean()
    print("5-fold CV AUC (PCA50, tumor vs normal):", aucs, "mean:", cv_auc_mean)

    names = list(rna_results.keys())
    pvals = [rna_results[n]["p_empirical"] for n in names]
    zscores = [rna_results[n]["z"] for n in names]
    rejected, p_adj = benjamini_hochberg(pvals, alpha=0.05)
    rows = [{"test": n, "real_max_H1": rna_results[n]["real"], "null_mean": rna_results[n]["null_mean"],
             "null_std": rna_results[n]["null_std"], "z_score": z, "p_empirical": p,
             "p_BH": padj, "reject_at_0.05_BH": rej, "passes_section7_threshold_z_gt_3": z is not None and z > 3.0}
            for n, z, p, padj, rej in zip(names, zscores, pvals, p_adj, rejected)]
    summary_df = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(summary_df.to_string(index=False))

    out_csv = os.path.join(_HERE, "tcga_luad_rnaseq_results_table.csv")
    summary_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # NOTE: real_max_h1_raw_hvg and real_max_h1_pca50 are computed entirely from the
    # real data (fixed seeds) and do not depend on n_gauss/n_perm.
    print("\n=== Verification against tcga_luad_report.md (RNA-seq layer) ===")
    checks = [
        ("real_max_h1_raw_hvg", real_raw_h1, TOL),
        ("real_max_h1_pca50", real_pca_h1, TOL),
        ("cv_auc_mean", cv_auc_mean, TOL_CV_AUC),
    ]
    for name, val, tol in checks:
        exp = EXPECTED[name]
        ok = abs(val - exp) < tol
        print(f"  {name}: got {val:.4f}, expected {exp:.4f}, tol={tol} -> {'PASS' if ok else 'FAIL'}")
    if not args.skip_assert:
        for name, val, tol in checks:
            exp = EXPECTED[name]
            assert abs(val - exp) < tol, f"{name} mismatch: got {val}, expected {exp} (tol {tol})"
        print("\nAll assertions PASSED.")
    if reduced:
        print("\nNOTE: reduced-draw-count run -- real-data numbers above are exact/assertable "
              "regardless, but null-distribution z-scores/p-values are noisier than the "
              "pre-registered n_gauss=500/n_perm=2000 run.")


if __name__ == "__main__":
    main()
