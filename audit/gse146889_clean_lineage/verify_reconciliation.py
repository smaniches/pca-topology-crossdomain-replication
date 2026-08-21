#!/usr/bin/env python3
"""Verify the matched-definition reconciliation for GSE146889.

This verifier is deliberately read-only. It independently reconstructs the
historical comparator quantities from committed clean-run machine-readable
results and asserts that the committed reconciled tables/claims use the same
definitions. It never changes raw, PCA, PH, null, bootstrap, or result files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

EXPECTED_INPUT_SHA256 = "16e22f1285cb3960cc30151e67f8bf1ccc94a5d2e42051de9cc47ca1bbba3fa4"
HISTORICAL = {
    "within_class": "results/confound_attribution_audit/table1_within_class_decomposition_GSE146889.csv",
    "within_stratum": "results/confound_attribution_audit/table2_within_stratum_control_GSE146889.csv",
    "residualization": "results/confound_attribution_audit/table3_residualization_control_GSE146889.csv",
    "bootstrap": "results/confound_attribution_audit/table4_block_bootstrap_ci_GSE146889.csv",
}


def close(a: float, b: float, tol: float = 1e-12) -> bool:
    return bool(np.isclose(float(a), float(b), rtol=0.0, atol=tol))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--results-dir", default="audit/gse146889_clean_lineage/results")
    args = parser.parse_args()
    repo = Path(args.repo_root).resolve()
    out = (repo / args.results_dir).resolve() if not Path(args.results_dir).is_absolute() else Path(args.results_dir)

    primary = json.loads((out / "primary_results.json").read_text(encoding="utf-8"))
    require(primary["data"]["sha256"] == EXPECTED_INPUT_SHA256, "unexpected GSE146889 input digest")
    require(primary["data"]["sample_count"] == 176, "unexpected sample count")
    require(primary["data"]["tumor_count"] == 91 and primary["data"]["normal_count"] == 85, "unexpected class counts")

    wc_new = pd.read_csv(out / "within_class.csv")
    wc_hist = pd.read_csv(repo / HISTORICAL["within_class"])
    for condition in sorted(set(wc_new.condition) & set(wc_hist.condition)):
        new = wc_new.loc[wc_new.condition == condition].iloc[0]
        hist = wc_hist.loc[wc_hist.condition == condition].iloc[0]
        require(close(new.observed_max_H1, hist.observed_max_H1), f"max-H1 mismatch for {condition}")

    ws_new = pd.read_csv(out / "within_stratum.csv")
    ws_hist = pd.read_csv(repo / HISTORICAL["within_stratum"])
    for quartile in ["Q1", "Q2", "Q3", "Q4"]:
        new = ws_new.loc[ws_new.quartile == quartile].iloc[0]
        hist = ws_hist.loc[ws_hist.quartile == quartile].iloc[0]
        require(close(new.observed_max_H1, hist.observed_max_H1), f"max-H1 mismatch for {quartile}")
        require(float(new.z_gauss) > 3.0 and float(new.z_perm) > 3.0, f"z>3 criterion failed for {quartile}")

    rz = pd.read_csv(out / "residualization.csv")
    hist_rz = pd.read_csv(repo / HISTORICAL["residualization"])
    real = primary["real_results"]["residualization"]
    mapping = {
        "Full HVG->PCA50 (class-mean intact)": (
            real["auc_C_1_intact"], real["auc_C_0_01_intact"], real["intact_cocycle_support"]
        ),
        "Class-mean-residualized HVG->PCA50": (
            real["auc_C_1_residualized"], real["auc_C_0_01_residualized"], real["residualized_cocycle_support"]
        ),
    }
    for space, (c1, c001, support) in mapping.items():
        row = rz.loc[rz.space == space].iloc[0]
        hist = hist_rz.loc[hist_rz.space == space].iloc[0]
        a = int(support["support_tumor"])
        b = int(support["support_normal"])
        c = 91 - a
        d = 85 - b
        p_two = float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
        p_greater = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        require(close(row.cv_auc_tumor_normal, c1["mean"]), f"C=1 AUC mismatch for {space}")
        require(close(row.cv_auc_tumor_normal, hist.cv_auc_tumor_normal), f"historical AUC mismatch for {space}")
        require(close(row.cv_auc_C_0_01_sensitivity, c001["mean"]), f"C=.01 sensitivity mismatch for {space}")
        require(close(row.fisher_p_tumor_enrichment, p_two), f"two-sided Fisher mismatch for {space}")
        require(close(row.fisher_p_tumor_enrichment, hist.fisher_p_tumor_enrichment), f"historical Fisher mismatch for {space}")
        require(close(row.fisher_p_tumor_enrichment_one_sided_sensitivity, p_greater), f"one-sided Fisher sensitivity mismatch for {space}")
        require(int(row.top_loop_n_participants) == int(support["support_n"]), f"cocycle support-n mismatch for {space}")
        require(int(row.top_loop_n_tumor) == a and int(row.top_loop_n_normal) == b, f"cocycle composition mismatch for {space}")

    residualized = rz.loc[rz.space == "Class-mean-residualized HVG->PCA50"].iloc[0]
    require(float(residualized.z_gaussian) > 3.0 and float(residualized.z_permutation) > 3.0, "residualized z criterion failed")

    bootstrap = pd.read_csv(out / "bootstrap.csv")
    mixed = bootstrap.loc[bootstrap.subset == "mixed (tumor+normal)"].iloc[0]
    tumor = bootstrap.loc[bootstrap.subset == "tumor-only"].iloc[0]
    normal = bootstrap.loc[bootstrap.subset == "normal-only"].iloc[0]
    require(bool(mixed.excludes_zero), "mixed bootstrap unexpectedly includes zero")
    require(not bool(tumor.excludes_zero) and not bool(normal.excludes_zero), "single-class bootstrap conclusion changed")

    claims = primary["claims"]
    require(bool(claims["deterministic_matched_definition_values_exact"]), "deterministic exact claim false")
    require(close(claims["deterministic_max_absolute_difference_vs_historical"], 0.0), "deterministic max difference is not zero")
    require(bool(claims["historical_qualitative_GSE146889_verdict_sustained"]), "qualitative verdict not sustained")
    require(bool(claims["aggregate_signal_survives_linear_residualization_z_gt_3_both_nulls"]), "aggregate residualization claim false")
    require(bool(claims["all_four_strata_exceed_z_3_against_both_nulls"]), "stratum robustness claim false")
    require(bool(claims["both_single_class_bootstrap_intervals_include_zero"]), "single-class bootstrap claim false")

    comparison = pd.read_csv(out / "historical_comparison.csv")
    deterministic = comparison.loc[comparison.metric_type == "deterministic"]
    require(len(deterministic) > 0, "no deterministic comparison rows")
    require((deterministic.status == "exact").all(), "a deterministic comparison row is not exact")
    require((pd.to_numeric(deterministic.absolute_difference, errors="coerce").fillna(np.inf) <= 1e-12).all(), "nonzero deterministic comparison difference")

    print(json.dumps({
        "input_sha256": EXPECTED_INPUT_SHA256,
        "samples": 176,
        "complete_matched_definition_deterministic_reproduction": True,
        "historical_qualitative_verdict_sustained": True,
        "reconciliation_verified_read_only": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
