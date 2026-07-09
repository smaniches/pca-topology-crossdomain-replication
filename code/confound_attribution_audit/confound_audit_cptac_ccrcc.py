"""
confound_audit_cptac_ccrcc.py -- confound-attribution audit for the CPTAC-CCRCC cohort
(clear-cell renal cell carcinoma proteomics, log2-ratio quantification), reproducing
results/confound_attribution_audit/confound_audit_CPTAC_CCRCC.md.

Data source: CPTAC Proteomic Data Commons (PDC) GraphQL API, study PDC000127
(log2_ratio quantification matrix + biospecimen metadata for sample_type).

Preprocessing (protein/proteomics-specific, differs from the RNA-seq cohorts):
  1. fillna(0.0) on the raw log2-ratio matrix
  2. median-center each SAMPLE (row) by its own median log2-ratio across all proteins
     (standard proteomics normalization, not applicable to RNA-seq FPKM/RPKM data)
  3. top-2000 most-variable proteins by variance across samples
  4. StandardScaler + PCA(50)

Topological backend: gudhi (matches what was actually run for this cohort in the
original audit -- gives identical H1 persistence to ripser for the same point cloud,
since both compute a Vietoris-Rips filtration, but gudhi was the one whose exact
numerical output the committed report's numbers come from).

Expected results (from the committed report): mixed max-H1 = 3.833, tumor-only max-H1 =
3.833 (EXACT match, same "tumor subset alone reproduces mixed magnitude" pattern as
GSE81089), residualized z_gaussian ~= 36.46.
"""
import os
import sys
import time
import json
import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
from confound_audit_common import (
    within_class_decomposition, confound_direction, within_stratum_control,
    residualization_control, block_bootstrap_ci,
)

PDC_URL = "https://pdc.cancer.gov/graphql"
PDC_STUDY_ID = "PDC000127"


def fetch_cptac_ccrcc(data_dir):
    """Fetch CPTAC-CCRCC log2-ratio proteomics matrix + sample-type labels from the PDC
    GraphQL API. Caches to data_dir so re-runs don't re-hit the API."""
    import requests
    os.makedirs(data_dir, exist_ok=True)
    raw_path = os.path.join(data_dir, "cptac_ccrcc_log2ratio_raw.csv")
    labels_path = os.path.join(data_dir, "cptac_ccrcc_labels.csv")
    if os.path.exists(raw_path) and os.path.exists(labels_path):
        return pd.read_csv(raw_path, index_col=0), pd.read_csv(labels_path, index_col=0)

    quant_query = ('{ quantDataMatrix( pdc_study_id: "' + PDC_STUDY_ID +
                    '" data_type: "log2_ratio" acceptDUA: true ) }')
    rq = requests.post(PDC_URL, json={"query": quant_query}, timeout=120)
    matrix = rq.json()["data"]["quantDataMatrix"]
    df = pd.DataFrame(matrix[1:], columns=matrix[0]).set_index("Gene/Aliquot")

    q_bio = """{ biospecimenPerStudy(pdc_study_id: "%s", acceptDUA: true) {
        aliquot_id sample_id case_id sample_submitter_id aliquot_submitter_id
        disease_type sample_type } }""" % PDC_STUDY_ID
    rb = requests.post(PDC_URL, json={"query": q_bio}, timeout=60)
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

    X_raw.to_csv(raw_path)
    labels.to_frame(name="sample_type").to_csv(labels_path)
    return X_raw, labels.to_frame(name="sample_type")


def preprocess(raw, labels):
    lab = labels.loc[raw.index, labels.columns[0]]
    tumor_mask = (lab == "Primary Tumor").values
    normal_mask = (lab == "Solid Tissue Normal").values

    X0 = raw.fillna(0.0)
    row_med = X0.median(axis=1)
    Xc = X0.sub(row_med, axis=0)
    var = Xc.var(axis=0, ddof=1)
    top2000 = var.sort_values(ascending=False).index[:2000]
    X_top = Xc[top2000].values
    scaler = StandardScaler()
    X_top_std = scaler.fit_transform(X_top)
    return X_top_std, tumor_mask, normal_mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.path.join(_ROOT, "data_cptac"))
    ap.add_argument("--results-dir", default=os.path.join(_ROOT, "results_cptac"))
    ap.add_argument("--n-gauss", type=int, default=300)
    ap.add_argument("--n-perm", type=int, default=300)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--quick", type=int, default=None)
    args = ap.parse_args()
    if args.quick:
        args.n_gauss = min(args.n_gauss, args.quick)
        args.n_perm = min(args.n_perm, args.quick)
        args.n_boot = min(args.n_boot, args.quick)
        print(f"[QUICK MODE] capping all draw counts to {args.quick}")

    os.makedirs(args.results_dir, exist_ok=True)
    t0 = time.time()
    raw, labels = fetch_cptac_ccrcc(args.data_dir)
    X_std, tumor_mask, normal_mask = preprocess(raw, labels)
    print(f"Loaded CPTAC-CCRCC: {X_std.shape}, tumor={tumor_mask.sum()}, normal={normal_mask.sum()} ({time.time()-t0:.0f}s)")

    BACKEND = "gudhi"

    # --- Control 1: within-class decomposition (raw standardized top-2000-protein space) ---
    print("\n=== Control 1: within-class decomposition ===")
    wc = within_class_decomposition(X_std, tumor_mask, n_gauss=args.n_gauss, n_perm=args.n_perm,
                                     n_features=2000, n_pca=None, backend=BACKEND, seed_base=1000)
    for name in ("mixed", "tumor", "normal"):
        r = wc[name]
        print(f"  {name:7s} n={r['n']:3d} observed={r['observed']:.4f} "
              f"z_gauss={r['gauss']['z']:.2f} z_perm={r['perm']['z']:.2f}")
    assert abs(wc["mixed"]["observed"] - wc["tumor"]["observed"]) < 1e-6, \
        "expected mixed and tumor-only max-H1 to match EXACTLY"
    print(f"  VERIFIED: mixed and tumor-only max-H1 match exactly ({wc['mixed']['observed']:.6f})")
    np.testing.assert_allclose(wc["mixed"]["observed"], 3.833, atol=5e-3)
    print("  VERIFIED: mixed/tumor-only max-H1 matches expected value 3.833 (tol 5e-3)")

    # --- Control 2: within-stratum quartile control (tumor samples, confound-projection bins) ---
    print("\n=== Control 2: within-stratum quartile control ===")
    u, proj_full = confound_direction(X_std, tumor_mask)
    tumor_idx = np.where(tumor_mask)[0]
    tumor_proj = proj_full[tumor_idx]
    q = pd.qcut(tumor_proj, 4, labels=False)
    X_by_stratum = [X_std[tumor_idx[q == i]] for i in range(4)]
    strata_df = within_stratum_control(X_by_stratum, n_gauss=min(args.n_gauss, 300), backend=BACKEND, seed_base=5000)
    print(strata_df[["quartile", "n", "observed", "gauss_z", "gauss_p"]].to_string(index=False))

    # --- Control 3: residualization ---
    print("\n=== Control 3: residualization ===")
    resid = residualization_control(X_std, tumor_mask, n_pcs=50, n_gauss=args.n_gauss,
                                     n_perm=args.n_perm, backend=BACKEND, seed_base=9000)
    print(f"  max abs class-mean diff after residualization: {resid['max_abs_classmean_diff_after_resid']:.2e}")
    print(f"  intact:       observed={resid['intact']['observed']:.4f} AUC={resid['intact']['auc']:.4f} "
          f"z_gauss={resid['intact']['gauss']['z']:.2f}")
    print(f"  residualized: observed={resid['residualized']['observed']:.4f} AUC={resid['residualized']['auc']:.4f} "
          f"z_gauss={resid['residualized']['gauss']['z']:.2f}")
    if not args.quick:
        np.testing.assert_allclose(resid['residualized']['gauss']['z'], 36.46, rtol=0.15)
        print("  VERIFIED: residualized z_gaussian matches expected ~36.46 (rtol 0.15)")
    else:
        print("  [QUICK MODE] skipping z_gaussian tolerance check (draw count reduced)")

    # --- Control 4: block-bootstrap CI (tumor-only) ---
    print("\n=== Control 4: block-bootstrap CI (tumor-only) ===")
    X_tumor = X_std[tumor_mask]
    null_mean = wc["tumor"]["gauss"]["null_mean"]
    null_sd = wc["tumor"]["gauss"]["null_sd"]
    boot = block_bootstrap_ci(X_tumor, null_mean, null_sd, n_boot=args.n_boot, seed=12345, backend=BACKEND)
    print(f"  delta 95% CI: [{boot['ci_delta'][0]:.2f}, {boot['ci_delta'][1]:.2f}] excludes_zero={boot['excludes_zero']}")
    print(f"  z 95% CI:     [{boot['ci_z'][0]:.2f}, {boot['ci_z'][1]:.2f}]")

    summary = pd.DataFrame([
        dict(control="within_class", subset="mixed", n=wc["mixed"]["n"], observed=wc["mixed"]["observed"], z_gauss=wc["mixed"]["gauss"]["z"]),
        dict(control="within_class", subset="tumor", n=wc["tumor"]["n"], observed=wc["tumor"]["observed"], z_gauss=wc["tumor"]["gauss"]["z"]),
        dict(control="within_class", subset="normal", n=wc["normal"]["n"], observed=wc["normal"]["observed"], z_gauss=wc["normal"]["gauss"]["z"]),
        dict(control="residualization", subset="intact", n=X_std.shape[0], observed=resid["intact"]["observed"], z_gauss=resid["intact"]["gauss"]["z"]),
        dict(control="residualization", subset="residualized", n=X_std.shape[0], observed=resid["residualized"]["observed"], z_gauss=resid["residualized"]["gauss"]["z"]),
    ])
    summary.to_csv(os.path.join(args.results_dir, "cptac_ccrcc_audit_summary.csv"), index=False)
    strata_df.to_csv(os.path.join(args.results_dir, "cptac_ccrcc_within_stratum.csv"), index=False)
    print(f"\nSaved results to {args.results_dir}")
    print(f"Total elapsed: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
