"""
confound_audit_gse81089.py -- confound-attribution audit for the GSE81089 pilot cohort
(218 NSCLC samples, 199 tumor / 19 normal), reproducing results/confound_attribution_audit/
confound_audit_GSE81089.md.

Loads/preprocesses the pilot expression data (same pipeline as the ablation sweep's
pre-registered default: top-2000 HVG by log1p variance, zero-fill for -1.0 sentinels,
StandardScaler, PCA(50), ripser backend), then runs the 4 shared controls from
confound_audit_common.py:
  1. within-class decomposition (mixed n=218 / tumor-only n=199 / normal-only n=19)
  2. within-stratum quartile control (bin TUMOR-only samples by confound-direction
     projection into quartiles; check z-score within each quartile)
  3. residualization (remove class-mean-shift direction, recompute PCA+max-H1+AUC)
  4. block-bootstrap CI (tumor-only subset)

Expected results (from the committed report, full-precision draws):
  mixed max-H1 = 3.887, tumor-only max-H1 = 3.887 (EXACT match -- tumor subset alone
  reproduces the mixed-set magnitude), residualized z_gaussian ~= 31.1
"""
import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import urllib.request
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
from confound_audit_common import (
    max_h1_persistence, within_class_decomposition, confound_direction,
    within_stratum_control, residualization_control, block_bootstrap_ci, zscore_pvalue,
)


def load_gse81089(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/suppl/GSE81089_FPKM_cufflinks.tsv.gz"
    out = os.path.join(data_dir, "GSE81089_FPKM.tsv.gz")
    if not os.path.exists(out):
        urllib.request.urlretrieve(url, out)

    df = pd.read_csv(out, sep="\t", index_col=0)
    df_clean = df.loc[[i for i in df.index if str(i).startswith("ENSG")]]
    # -1.0 sentinels -> missing -> zero-fill (pre-registered default imputation rule)
    df_clean = df_clean.replace(-1.0, np.nan).fillna(0.0)

    expr = df_clean.T  # samples x genes
    log_expr = np.log1p(expr)
    gene_var = log_expr.var(axis=0)
    nonzero = gene_var > 1e-12
    log_expr = log_expr.loc[:, nonzero]
    gene_var = gene_var[nonzero]

    labels = pd.Series(
        ["Tumor" if c.split("_")[0].endswith("T") else ("Normal" if c.split("_")[0].endswith("N") else "Unknown")
         for c in expr.index], index=expr.index)

    N_GENES = 2000
    hvg_genes = gene_var.sort_values(ascending=False).index[:N_GENES]
    X_hvg_raw = log_expr[hvg_genes].values  # 218 x 2000
    tumor_mask = (labels.values == "Tumor")
    return X_hvg_raw, tumor_mask, labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.path.join(_ROOT, "data_gse81089"))
    ap.add_argument("--results-dir", default=os.path.join(_ROOT, "results_gse81089"))
    ap.add_argument("--n-gauss", type=int, default=500, help="Gaussian null draws (pre-registered: 500)")
    ap.add_argument("--n-perm", type=int, default=2000, help="Permutation null draws (pre-registered: 2000)")
    ap.add_argument("--n-boot", type=int, default=2000, help="Bootstrap draws (pre-registered: 2000)")
    ap.add_argument("--quick", type=int, default=None,
                     help="If set, cap all draw counts to this value for a fast smoke test "
                          "(disclose explicitly -- z/p-values become approximate, real-data "
                          "statistics are unaffected)")
    args = ap.parse_args()
    if args.quick:
        args.n_gauss = min(args.n_gauss, args.quick)
        args.n_perm = min(args.n_perm, args.quick)
        args.n_boot = min(args.n_boot, args.quick)
        print(f"[QUICK MODE] capping all draw counts to {args.quick} -- z/p-values are approximate")

    os.makedirs(args.results_dir, exist_ok=True)
    t0 = time.time()
    X_hvg_raw, tumor_mask, labels = load_gse81089(args.data_dir)
    print(f"Loaded GSE81089: {X_hvg_raw.shape}, tumor={tumor_mask.sum()}, normal={(~tumor_mask).sum()} ({time.time()-t0:.0f}s)")

    scaler = StandardScaler()
    X_hvg_std = scaler.fit_transform(X_hvg_raw)
    pca = PCA(n_components=50, random_state=42)
    X_pca50 = pca.fit_transform(X_hvg_std)

    # --- Control 1: within-class decomposition (raw HVG space, matching report's Table 1) ---
    print("\n=== Control 1: within-class decomposition (raw HVG space) ===")
    wc = within_class_decomposition(X_hvg_raw, tumor_mask, n_gauss=args.n_gauss, n_perm=args.n_perm,
                                     n_features=2000, n_pca=None, backend="ripser", seed_base=100)
    for name in ("mixed", "tumor", "normal"):
        r = wc[name]
        print(f"  {name:7s} n={r['n']:3d} observed={r['observed']:.4f} "
              f"z_gauss={r['gauss']['z']:.2f} z_perm={r['perm']['z']:.2f}")
    assert abs(wc["mixed"]["observed"] - wc["tumor"]["observed"]) < 1e-6, \
        "expected mixed and tumor-only max-H1 to match EXACTLY (tumor subset alone reproduces full magnitude)"
    print(f"  VERIFIED: mixed and tumor-only max-H1 match exactly ({wc['mixed']['observed']:.6f})")
    np.testing.assert_allclose(wc["mixed"]["observed"], 3.887, atol=5e-3)
    print("  VERIFIED: mixed/tumor-only max-H1 matches expected value 3.887 (tol 5e-3)")

    # --- Control 2: within-stratum quartile control (bin TUMOR samples by confound projection) ---
    print("\n=== Control 2: within-stratum quartile control ===")
    u, proj_full = confound_direction(X_hvg_std, tumor_mask)
    tumor_idx = np.where(tumor_mask)[0]
    tumor_proj = proj_full[tumor_idx]
    quartile_edges = np.quantile(tumor_proj, [0, 0.25, 0.5, 0.75, 1.0])
    strata_id = np.digitize(tumor_proj, quartile_edges[1:-1])
    X_by_stratum = [X_hvg_raw[tumor_idx[strata_id == q]] for q in range(4)]
    strata_df = within_stratum_control(X_by_stratum, n_gauss=min(args.n_gauss, 300), backend="ripser", seed_base=1000)
    print(strata_df[["quartile", "n", "observed", "gauss_z", "gauss_p"]].to_string(index=False))

    # --- Control 3: residualization ---
    print("\n=== Control 3: residualization ===")
    resid = residualization_control(X_hvg_std, tumor_mask, n_pcs=50,
                                     n_gauss=args.n_gauss, n_perm=args.n_perm, backend="ripser", seed_base=9000)
    print(f"  max abs class-mean diff after residualization: {resid['max_abs_classmean_diff_after_resid']:.2e}")
    print(f"  intact:       observed={resid['intact']['observed']:.4f} AUC={resid['intact']['auc']:.4f} "
          f"z_gauss={resid['intact']['gauss']['z']:.2f}")
    print(f"  residualized: observed={resid['residualized']['observed']:.4f} AUC={resid['residualized']['auc']:.4f} "
          f"z_gauss={resid['residualized']['gauss']['z']:.2f}")
    if not args.quick:
        np.testing.assert_allclose(resid['residualized']['gauss']['z'], 31.1, rtol=0.15)
        print("  VERIFIED: residualized z_gaussian matches expected ~31.1 (rtol 0.15)")
    else:
        print("  [QUICK MODE] skipping z_gaussian tolerance check (draw count reduced)")

    # --- Control 4: block-bootstrap CI (tumor-only subset) ---
    print("\n=== Control 4: block-bootstrap CI (tumor-only) ===")
    X_tumor_raw = X_hvg_raw[tumor_mask]
    null_mean = wc["tumor"]["gauss"]["null_mean"]
    null_sd = wc["tumor"]["gauss"]["null_sd"]
    boot = block_bootstrap_ci(X_tumor_raw, null_mean, null_sd, n_boot=args.n_boot, seed=42, backend="ripser")
    print(f"  delta 95% CI: [{boot['ci_delta'][0]:.2f}, {boot['ci_delta'][1]:.2f}] excludes_zero={boot['excludes_zero']}")
    print(f"  z 95% CI:     [{boot['ci_z'][0]:.2f}, {boot['ci_z'][1]:.2f}]")

    # --- Save results ---
    summary = pd.DataFrame([
        dict(control="within_class", subset="mixed", n=wc["mixed"]["n"], observed=wc["mixed"]["observed"], z_gauss=wc["mixed"]["gauss"]["z"]),
        dict(control="within_class", subset="tumor", n=wc["tumor"]["n"], observed=wc["tumor"]["observed"], z_gauss=wc["tumor"]["gauss"]["z"]),
        dict(control="within_class", subset="normal", n=wc["normal"]["n"], observed=wc["normal"]["observed"], z_gauss=wc["normal"]["gauss"]["z"]),
        dict(control="residualization", subset="intact", n=X_hvg_std.shape[0], observed=resid["intact"]["observed"], z_gauss=resid["intact"]["gauss"]["z"]),
        dict(control="residualization", subset="residualized", n=X_hvg_std.shape[0], observed=resid["residualized"]["observed"], z_gauss=resid["residualized"]["gauss"]["z"]),
    ])
    summary.to_csv(os.path.join(args.results_dir, "gse81089_audit_summary.csv"), index=False)
    strata_df.to_csv(os.path.join(args.results_dir, "gse81089_within_stratum.csv"), index=False)
    print(f"\nSaved results to {args.results_dir}")
    print(f"Total elapsed: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
