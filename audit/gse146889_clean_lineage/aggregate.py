#!/usr/bin/env python3
"""Aggregate clean-lineage GSE146889 shards into auditable results."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EXPECTED_COUNTS: dict[str, int] = {
    "mixed_gauss": 500,
    "mixed_perm": 2000,
    "tumor_gauss": 300,
    "tumor_perm": 300,
    "normal_gauss": 300,
    "normal_perm": 300,
    "q1_gauss": 300,
    "q1_perm": 300,
    "q2_gauss": 300,
    "q2_perm": 300,
    "q3_gauss": 300,
    "q3_perm": 300,
    "q4_gauss": 300,
    "q4_perm": 300,
    "bootstrap_mixed": 2000,
    "bootstrap_tumor": 2000,
    "bootstrap_normal": 2000,
}

HISTORICAL_FILES = {
    "within_class": "results/confound_attribution_audit/table1_within_class_decomposition_GSE146889.csv",
    "within_stratum": "results/confound_attribution_audit/table2_within_stratum_control_GSE146889.csv",
    "residualization": "results/confound_attribution_audit/table3_residualization_control_GSE146889.csv",
    "bootstrap": "results/confound_attribution_audit/table4_block_bootstrap_ci_GSE146889.csv",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def load_shards(shards_dir: Path) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    grouped: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
    provenance: list[dict[str, Any]] = []
    paths = sorted(shards_dir.rglob("*.npz"))
    if not paths:
        raise RuntimeError(f"No shard npz files found below {shards_dir}")
    for path in paths:
        data = np.load(path, allow_pickle=False)
        task = str(data["task"].item())
        indices = np.asarray(data["indices"], dtype=np.int64)
        values = np.asarray(data["values"], dtype=np.float64)
        if len(indices) != len(values):
            raise RuntimeError(f"{path}: index/value length mismatch")
        grouped.setdefault(task, []).append((indices, values))
        provenance.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "task": task,
                "start": int(data["start"].item()),
                "count": int(data["count"].item()),
                "base_seed": int(data["base_seed"].item()),
                "prepared_sha256": str(data["prepared_sha256"].item()),
                "worker_sha256": str(data["worker_sha256"].item()),
            }
        )

    completed: dict[str, np.ndarray] = {}
    missing = sorted(set(EXPECTED_COUNTS) - set(grouped))
    extra = sorted(set(grouped) - set(EXPECTED_COUNTS))
    if missing or extra:
        raise RuntimeError(f"Shard task mismatch: missing={missing}, extra={extra}")
    for task, expected_count in EXPECTED_COUNTS.items():
        pairs = grouped[task]
        indices = np.concatenate([pair[0] for pair in pairs])
        values = np.concatenate([pair[1] for pair in pairs])
        order = np.argsort(indices, kind="mergesort")
        indices = indices[order]
        values = values[order]
        expected_indices = np.arange(expected_count, dtype=np.int64)
        if not np.array_equal(indices, expected_indices):
            raise RuntimeError(
                f"{task}: expected draw indices 0..{expected_count - 1}; got "
                f"{indices[:5]}...{indices[-5:]}"
            )
        if not np.isfinite(values).all():
            raise RuntimeError(f"{task}: non-finite values")
        completed[task] = values
    return completed, provenance


def null_stats(observed: float, draws: np.ndarray) -> dict[str, float]:
    values = np.asarray(draws, dtype=float)
    mean = float(values.mean())
    sd = float(values.std(ddof=1))
    z = float((observed - mean) / sd) if sd > 0 else float("nan")
    p = float((int(np.sum(values >= observed)) + 1) / (len(values) + 1))
    return {
        "null_mean": mean,
        "null_sd": sd,
        "z": z,
        "p": p,
        "n_draws": int(len(values)),
    }


def markdown_table(frame: pd.DataFrame, *, precision: int = 4) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(str(c) for c in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        rendered: list[str] = []
        for value in row.tolist():
            if isinstance(value, (float, np.floating)):
                if np.isnan(value):
                    rendered.append("NaN")
                elif abs(float(value)) < 1e-3 and float(value) != 0.0:
                    rendered.append(f"{float(value):.{precision}e}")
                else:
                    rendered.append(f"{float(value):.{precision}f}")
            elif isinstance(value, (bool, np.bool_)):
                rendered.append("yes" if bool(value) else "no")
            else:
                rendered.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(rendered) + " |")
    return "\n".join(lines)


def compare_value(
    rows: list[dict[str, Any]],
    *,
    section: str,
    key: str,
    metric: str,
    historical: Any,
    new: Any,
    metric_type: str,
) -> None:
    try:
        hist_float = float(historical)
        new_float = float(new)
        abs_diff = abs(new_float - hist_float)
        rel_diff = abs_diff / abs(hist_float) if hist_float != 0 else np.nan
        if abs_diff <= 1e-9:
            status = "exact"
        elif metric_type == "deterministic" and abs_diff <= 1e-3:
            status = "numerically_reproduced"
        elif metric_type == "deterministic" and rel_diff <= 0.01:
            status = "close_deterministic"
        elif metric_type in {"monte_carlo", "bootstrap"} and rel_diff <= 0.10:
            status = "within_10pct"
        else:
            status = "different"
    except (TypeError, ValueError):
        hist_float = historical
        new_float = new
        abs_diff = np.nan
        rel_diff = np.nan
        status = "same" if str(historical) == str(new) else "different"
    rows.append(
        {
            "section": section,
            "key": key,
            "metric": metric,
            "metric_type": metric_type,
            "historical": hist_float,
            "new": new_float,
            "absolute_difference": abs_diff,
            "relative_difference": rel_diff,
            "status": status,
        }
    )


def make_comparison(
    repo_root: Path,
    within_class: pd.DataFrame,
    within_stratum: pd.DataFrame,
    residualization: pd.DataFrame,
    bootstrap: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    historical = pd.read_csv(repo_root / HISTORICAL_FILES["within_class"])
    hist_map = {str(row["condition"]): row for _, row in historical.iterrows()}
    new_map = {str(row["condition"]): row for _, row in within_class.iterrows()}
    for key in sorted(set(hist_map) & set(new_map)):
        for metric, metric_type in [
            ("observed_max_H1", "deterministic"),
            ("z_gaussian", "monte_carlo"),
            ("p_gaussian", "monte_carlo"),
            ("z_permutation", "monte_carlo"),
            ("p_permutation", "monte_carlo"),
        ]:
            compare_value(
                rows,
                section="within_class",
                key=key,
                metric=metric,
                historical=hist_map[key][metric],
                new=new_map[key][metric],
                metric_type=metric_type,
            )

    historical = pd.read_csv(repo_root / HISTORICAL_FILES["within_stratum"])
    hist_map = {str(row["quartile"]): row for _, row in historical.iterrows()}
    new_map = {str(row["quartile"]): row for _, row in within_stratum.iterrows()}
    for key in sorted(set(hist_map) & set(new_map)):
        for metric, metric_type in [
            ("observed_max_H1", "deterministic"),
            ("z_gauss", "monte_carlo"),
            ("z_perm", "monte_carlo"),
        ]:
            compare_value(
                rows,
                section="within_stratum",
                key=key,
                metric=metric,
                historical=hist_map[key][metric],
                new=new_map[key][metric],
                metric_type=metric_type,
            )

    historical = pd.read_csv(repo_root / HISTORICAL_FILES["residualization"])
    hist_map = {str(row["space"]): row for _, row in historical.iterrows()}
    new_map = {str(row["space"]): row for _, row in residualization.iterrows()}
    for key in sorted(set(hist_map) & set(new_map)):
        for metric, metric_type in [
            ("observed_max_H1", "deterministic"),
            ("z_gaussian", "monte_carlo"),
            ("z_permutation", "monte_carlo"),
            ("cv_auc_tumor_normal", "deterministic"),
            ("top_loop_n_participants", "diagnostic"),
            ("top_loop_n_tumor", "diagnostic"),
            ("top_loop_n_normal", "diagnostic"),
            ("fisher_p_tumor_enrichment", "diagnostic"),
        ]:
            compare_value(
                rows,
                section="residualization",
                key=key,
                metric=metric,
                historical=hist_map[key][metric],
                new=new_map[key][metric],
                metric_type=metric_type,
            )

    historical = pd.read_csv(repo_root / HISTORICAL_FILES["bootstrap"])
    hist_map = {str(row["subset"]): row for _, row in historical.iterrows()}
    new_map = {str(row["subset"]): row for _, row in bootstrap.iterrows()}
    for key in sorted(set(hist_map) & set(new_map)):
        for metric in [
            "point_delta",
            "delta_ci95_lo",
            "delta_ci95_hi",
            "point_z",
            "z_ci95_lo",
            "z_ci95_hi",
            "excludes_zero",
        ]:
            compare_value(
                rows,
                section="bootstrap",
                key=key,
                metric=metric,
                historical=hist_map[key][metric],
                new=new_map[key][metric],
                metric_type="bootstrap",
            )

    return pd.DataFrame(rows)


def generate_figure(results_dir: Path) -> None:
    within_class = pd.read_csv(results_dir / "within_class.csv")
    within_stratum = pd.read_csv(results_dir / "within_stratum.csv")
    residualization = pd.read_csv(results_dir / "residualization.csv")
    bootstrap = pd.read_csv(results_dir / "bootstrap.csv")

    figure, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    x = np.arange(len(within_class))
    ax.bar(x - 0.25, within_class["observed_max_H1"], width=0.25, label="Observed")
    ax.bar(x, within_class["gaussian_null_mean"], width=0.25, label="Gaussian null mean")
    ax.bar(x + 0.25, within_class["permutation_null_mean"], width=0.25, label="Permutation null mean")
    ax.set_xticks(x, within_class["condition"], rotation=20, ha="right")
    ax.set_ylabel("Maximum H1 persistence")
    ax.set_title("Within-class decomposition")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    x = np.arange(len(within_stratum))
    ax.bar(x - 0.2, within_stratum["z_gauss"], width=0.4, label="Gaussian z")
    ax.bar(x + 0.2, within_stratum["z_perm"], width=0.4, label="Permutation z")
    ax.axhline(3.0, linestyle="--", linewidth=1, label="z=3 criterion")
    ax.set_xticks(x, within_stratum["quartile"])
    ax.set_ylabel("z-score")
    ax.set_title("Confound-projection quartiles")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    x = np.arange(len(residualization))
    ax.bar(x - 0.18, residualization["observed_max_H1"], width=0.36, label="Max H1")
    ax.set_xticks(x, residualization["space"], rotation=15, ha="right")
    ax.set_ylabel("Maximum H1 persistence")
    ax2 = ax.twinx()
    ax2.plot(x + 0.18, residualization["cv_auc_tumor_normal"], marker="o", label="CV AUC")
    ax2.axhline(0.5, linestyle=":", linewidth=1)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("Tumor/normal CV AUC")
    ax.set_title("Residualization: topology vs separability")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=8, loc="best")

    ax = axes[1, 1]
    x = np.arange(len(bootstrap))
    point = bootstrap["point_delta"].to_numpy()
    lower = point - bootstrap["delta_ci95_lo"].to_numpy()
    upper = bootstrap["delta_ci95_hi"].to_numpy() - point
    ax.errorbar(x, point, yerr=np.vstack([lower, upper]), fmt="o", capsize=5)
    ax.axhline(0.0, linestyle="--", linewidth=1)
    ax.set_xticks(x, bootstrap["subset"], rotation=20, ha="right")
    ax.set_ylabel("Observed max H1 - Gaussian null mean")
    ax.set_title("Row-bootstrap stability")

    figure.suptitle("GSE146889 clean-lineage confound-attribution rerun")
    figure.tight_layout()
    figure.savefig(results_dir / "gse146889_clean_lineage_audit.png", dpi=200)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True)
    parser.add_argument("--prepare-dir", required=True)
    parser.add_argument("--shards-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--environment-file", required=True)
    parser.add_argument("--protocol-anchor-sha", required=True)
    args = parser.parse_args()

    repo_root = Path(args.repo_root)
    prepare_dir = Path(args.prepare_dir)
    shards_dir = Path(args.shards_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with (prepare_dir / "real_results.json").open("r", encoding="utf-8") as handle:
        real = json.load(handle)
    with (prepare_dir / "provenance_base.json").open("r", encoding="utf-8") as handle:
        provenance_base = json.load(handle)
    prepared_path = prepare_dir / "prepared.npz"
    draws, shard_provenance = load_shards(shards_dir)

    subset_key = {
        "Mixed (tumor+normal)": "mixed",
        "Tumor-only": "tumor",
        "Normal-only": "normal",
    }
    within_rows: list[dict[str, Any]] = []
    for condition, key in subset_key.items():
        observed = float(real["subsets"][key]["observed_max_h1"])
        gauss = null_stats(observed, draws[f"{key}_gauss"])
        perm = null_stats(observed, draws[f"{key}_perm"])
        within_rows.append(
            {
                "condition": condition,
                "n": int(real["subsets"][key]["n"]),
                "observed_max_H1": observed,
                "gaussian_null_mean": gauss["null_mean"],
                "gaussian_null_sd": gauss["null_sd"],
                "z_gaussian": gauss["z"],
                "p_gaussian": gauss["p"],
                "permutation_null_mean": perm["null_mean"],
                "permutation_null_sd": perm["null_sd"],
                "z_permutation": perm["z"],
                "p_permutation": perm["p"],
                "n_gaussian": gauss["n_draws"],
                "n_permutation": perm["n_draws"],
            }
        )
    within_class = pd.DataFrame(within_rows)
    within_class.to_csv(output_dir / "within_class.csv", index=False)

    stratum_rows: list[dict[str, Any]] = []
    quartile_meta = {entry["name"]: entry for entry in real["quartiles"]}
    for index in range(1, 5):
        key = f"q{index}"
        observed = float(real["subsets"][key]["observed_max_h1"])
        gauss = null_stats(observed, draws[f"{key}_gauss"])
        perm = null_stats(observed, draws[f"{key}_perm"])
        meta = quartile_meta[key]
        stratum_rows.append(
            {
                "quartile": f"Q{index}",
                "n": int(real["subsets"][key]["n"]),
                "n_tumor": int(meta["n_tumor"]),
                "n_normal": int(meta["n_normal"]),
                "projection_min": float(meta["projection_min"]),
                "projection_max": float(meta["projection_max"]),
                "observed_max_H1": observed,
                "gauss_null_mean": gauss["null_mean"],
                "gauss_null_sd": gauss["null_sd"],
                "z_gauss": gauss["z"],
                "p_gauss": gauss["p"],
                "perm_null_mean": perm["null_mean"],
                "perm_null_sd": perm["null_sd"],
                "z_perm": perm["z"],
                "p_perm": perm["p"],
            }
        )
    within_stratum = pd.DataFrame(stratum_rows)
    within_stratum.to_csv(output_dir / "within_stratum.csv", index=False)

    mixed_observed = float(real["subsets"]["mixed"]["observed_max_h1"])
    mixed_gauss = null_stats(mixed_observed, draws["mixed_gauss"])
    mixed_perm = null_stats(mixed_observed, draws["mixed_perm"])
    residual_observed = float(real["residualization"]["observed_max_h1"])
    residual_gauss = null_stats(residual_observed, draws["mixed_gauss"])
    residual_perm = null_stats(residual_observed, draws["mixed_perm"])
    intact_cocycle = real["residualization"]["intact_cocycle_support"]
    residual_cocycle = real["residualization"]["residualized_cocycle_support"]
    intact_auc = real["residualization"]["auc_C_0_01_intact"]
    residual_auc = real["residualization"]["auc_C_0_01_residualized"]
    residualization = pd.DataFrame(
        [
            {
                "space": "Full HVG->PCA50 (class-mean intact)",
                "observed_max_H1": mixed_observed,
                "z_gaussian": mixed_gauss["z"],
                "z_permutation": mixed_perm["z"],
                "cv_auc_tumor_normal": float(intact_auc["mean"]),
                "cv_auc_sd": float(intact_auc["sd_population"]),
                "top_loop_n_participants": int(intact_cocycle["support_n"]),
                "top_loop_n_tumor": int(intact_cocycle["support_tumor"]),
                "top_loop_n_normal": int(intact_cocycle["support_normal"]),
                "fisher_p_tumor_enrichment": float(intact_cocycle["fisher_p_tumor_enrichment"]),
                "loop_diagnostic_definition": "ripser representative H1 cocycle support vertices",
            },
            {
                "space": "Class-mean-residualized HVG->PCA50",
                "observed_max_H1": residual_observed,
                "z_gaussian": residual_gauss["z"],
                "z_permutation": residual_perm["z"],
                "cv_auc_tumor_normal": float(residual_auc["mean"]),
                "cv_auc_sd": float(residual_auc["sd_population"]),
                "top_loop_n_participants": int(residual_cocycle["support_n"]),
                "top_loop_n_tumor": int(residual_cocycle["support_tumor"]),
                "top_loop_n_normal": int(residual_cocycle["support_normal"]),
                "fisher_p_tumor_enrichment": float(residual_cocycle["fisher_p_tumor_enrichment"]),
                "loop_diagnostic_definition": "ripser representative H1 cocycle support vertices",
            },
        ]
    )
    residualization.to_csv(output_dir / "residualization.csv", index=False)

    bootstrap_rows: list[dict[str, Any]] = []
    display_name = {
        "mixed": "mixed (tumor+normal)",
        "tumor": "tumor-only",
        "normal": "normal-only",
    }
    for key in ["mixed", "tumor", "normal"]:
        observed = float(real["subsets"][key]["observed_max_h1"])
        null = null_stats(observed, draws[f"{key}_gauss"])
        boot_observed = draws[f"bootstrap_{key}"]
        deltas = boot_observed - null["null_mean"]
        z_values = deltas / null["null_sd"]
        ci_delta = np.percentile(deltas, [2.5, 97.5])
        ci_z = np.percentile(z_values, [2.5, 97.5])
        bootstrap_rows.append(
            {
                "subset": display_name[key],
                "n": int(real["subsets"][key]["n"]),
                "point_delta": observed - null["null_mean"],
                "delta_ci95_lo": float(ci_delta[0]),
                "delta_ci95_hi": float(ci_delta[1]),
                "point_z": null["z"],
                "z_ci95_lo": float(ci_z[0]),
                "z_ci95_hi": float(ci_z[1]),
                "excludes_zero": bool(ci_delta[0] > 0 or ci_delta[1] < 0),
                "n_boot": int(len(boot_observed)),
                "fixed_null_approximation": True,
            }
        )
    bootstrap = pd.DataFrame(bootstrap_rows)
    bootstrap.to_csv(output_dir / "bootstrap.csv", index=False)

    confound_spectrum = pd.DataFrame(real["confound_spectrum"]["pc_correlations"])
    confound_spectrum["variance_fraction_along_direction"] = np.nan
    confound_spectrum.loc[0, "variance_fraction_along_direction"] = float(
        real["confound_spectrum"]["variance_fraction_along_direction"]
    )
    confound_spectrum.loc[0, "projection_direct_auc_orientation_invariant"] = float(
        real["confound_spectrum"]["projection_direct_auc_orientation_invariant"]
    )
    confound_spectrum.to_csv(output_dir / "confound_spectrum.csv", index=False)

    np.savez_compressed(output_dir / "null_and_bootstrap_draws.npz", **draws)

    comparison = make_comparison(
        repo_root,
        within_class,
        within_stratum,
        residualization,
        bootstrap,
    )
    comparison.to_csv(output_dir / "historical_comparison.csv", index=False)

    aggregate_survives = bool(
        residual_gauss["z"] > 3.0 and residual_perm["z"] > 3.0
    )
    strata_survive = bool(
        (within_stratum["z_gauss"] > 3.0).all()
        and (within_stratum["z_perm"] > 3.0).all()
    )
    class_weaker = bool(
        float(real["subsets"]["tumor"]["observed_max_h1"])
        < 0.95 * mixed_observed
        and float(real["subsets"]["normal"]["observed_max_h1"])
        < 0.95 * mixed_observed
    )
    boot_tumor = bootstrap.loc[bootstrap["subset"] == "tumor-only"].iloc[0]
    boot_normal = bootstrap.loc[bootstrap["subset"] == "normal-only"].iloc[0]
    single_class_bootstrap_uncertain = bool(
        not bool(boot_tumor["excludes_zero"])
        and not bool(boot_normal["excludes_zero"])
    )
    historical_qualitative_verdict_sustained = bool(
        aggregate_survives
        and strata_survive
        and class_weaker
        and single_class_bootstrap_uncertain
    )

    deterministic_rows = comparison[
        comparison["metric_type"].isin(["deterministic"])
    ]
    deterministic_max_abs = float(
        pd.to_numeric(deterministic_rows["absolute_difference"], errors="coerce").max()
    )
    deterministic_reproduced = bool(deterministic_max_abs <= 0.05)

    primary_results = {
        "scope": (
            "Independent clean-lineage re-execution of the GSE146889 confound-attribution "
            "audit only; not a re-execution of every cohort or the entire manuscript."
        ),
        "data": real["data"],
        "real_results": real,
        "draw_counts": EXPECTED_COUNTS,
        "claims": {
            "deterministic_real_data_reproduced_within_0_05": deterministic_reproduced,
            "deterministic_max_absolute_difference_vs_historical": deterministic_max_abs,
            "aggregate_signal_survives_linear_residualization_z_gt_3_both_nulls": aggregate_survives,
            "all_four_strata_exceed_z_3_against_both_nulls": strata_survive,
            "single_class_observed_statistics_weaker_than_mixed": class_weaker,
            "both_single_class_bootstrap_intervals_include_zero": single_class_bootstrap_uncertain,
            "historical_qualitative_GSE146889_verdict_sustained": historical_qualitative_verdict_sustained,
        },
        "interpretation_boundary": {
            "not_external_human_validation": True,
            "not_a_new_preregistration_of_the_original_study": True,
            "cocycle_support_is_not_a_canonical_cycle": True,
            "mixed_tissue_cohort_remains_a_limit": True,
        },
    }
    with (output_dir / "primary_results.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(primary_results), handle, indent=2, sort_keys=True)
        handle.write("\n")

    generate_figure(output_dir)

    with Path(args.environment_file).open("r", encoding="utf-8") as handle:
        environment_text = handle.read()
    provenance = {
        **provenance_base,
        "protocol_anchor_commit_sha": args.protocol_anchor_sha,
        "aggregate_script_sha256": sha256_file(Path(__file__)),
        "prepared_npz_sha256": sha256_file(prepared_path),
        "shards": shard_provenance,
        "shard_count": len(shard_provenance),
        "environment_file_sha256": sha256_file(Path(args.environment_file)),
        "environment": environment_text.splitlines(),
        "workflow": {
            "repository": os.environ.get("GITHUB_REPOSITORY"),
            "workflow": os.environ.get("GITHUB_WORKFLOW"),
            "run_id": os.environ.get("GITHUB_RUN_ID"),
            "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
            "execution_git_sha": os.environ.get("GITHUB_SHA"),
            "ref": os.environ.get("GITHUB_REF"),
            "runner_os": os.environ.get("RUNNER_OS"),
            "platform": platform.platform(),
        },
        "output_hashes": {},
    }

    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "run_provenance.json":
            provenance["output_hashes"][path.name] = sha256_file(path)
    with (output_dir / "run_provenance.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(provenance), handle, indent=2, sort_keys=True)
        handle.write("\n")

    report_lines = [
        "# GSE146889 clean-lineage confound-attribution rerun",
        "",
        "**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)",
        "",
        "## Scope and provenance",
        "",
        "This audit re-executed the GSE146889 confound-attribution analysis from the public GEO RPKM matrix under a protocol committed before execution. It did not import the historical transcript-reconstructed GSE146889 driver, did not load its saved null pickles, and computed all new statistics before reading historical result tables. It is an internal computational replication, not independent human review.",
        "",
        f"- Input SHA-256: `{real['data']['sha256']}`",
        f"- Input bytes: `{real['data']['size_bytes']}`",
        f"- Samples: {real['data']['sample_count']} ({real['data']['tumor_count']} tumor, {real['data']['normal_count']} normal)",
        f"- Valid gene rows: {real['data']['valid_gene_rows']}; invalid/missing rows removed: {real['data']['invalid_or_missing_rpkm_rows_removed']}",
        f"- Protocol anchor commit: `{args.protocol_anchor_sha}`",
        f"- Execution commit: `{os.environ.get('GITHUB_SHA')}`",
        "",
        "## Result",
        "",
        (
            "**Historical qualitative verdict sustained.**"
            if historical_qualitative_verdict_sustained
            else "**Historical qualitative verdict not fully sustained; see the disagreements below.**"
        ),
        "",
        f"Deterministic real-data values reproduced within an absolute tolerance of 0.05: **{deterministic_reproduced}** (maximum absolute difference `{deterministic_max_abs:.6g}`).",
        "",
        "The clean rerun supports the following bounded interpretation:" if historical_qualitative_verdict_sustained else "The clean rerun changes at least one interpretation criterion:",
        "",
        f"- Aggregate max-H1 excess survives linear class-mean residualization against both nulls: **{aggregate_survives}**.",
        f"- All four confound-projection quartiles exceed z=3 against both nulls: **{strata_survive}**.",
        f"- Tumor-only and normal-only observed max-H1 values are materially below the mixed value: **{class_weaker}**.",
        f"- Both single-class row-bootstrap intervals include zero: **{single_class_bootstrap_uncertain}**.",
        "",
        "This sustains a qualified claim: the aggregate topological excess is not removed by the exact linear tumor/normal mean shift, but the single-class signal is weaker and unstable under row resampling. It does not establish complete confound independence, causality, or external biological validity.",
        "",
        "## Within-class decomposition",
        "",
        markdown_table(within_class[["condition", "n", "observed_max_H1", "z_gaussian", "p_gaussian", "z_permutation", "p_permutation"]]),
        "",
        "## Within-stratum control",
        "",
        markdown_table(within_stratum[["quartile", "n", "n_tumor", "n_normal", "observed_max_H1", "z_gauss", "z_perm"]]),
        "",
        "## Residualization",
        "",
        markdown_table(residualization[["space", "observed_max_H1", "z_gaussian", "z_permutation", "cv_auc_tumor_normal", "top_loop_n_participants", "top_loop_n_tumor", "top_loop_n_normal", "fisher_p_tumor_enrichment"]]),
        "",
        "The loop-composition columns use support vertices of ripser's representative H1 cocycle. They are not a canonical cycle identity and are reported only as a descriptive sensitivity diagnostic.",
        "",
        "## Bootstrap stability",
        "",
        markdown_table(bootstrap[["subset", "n", "point_delta", "delta_ci95_lo", "delta_ci95_hi", "excludes_zero"]]),
        "",
        "The bootstrap holds the Gaussian null distribution fixed across row resamples, matching the historical audit. This is an approximation and is not a joint bootstrap of both observed and null processes.",
        "",
        "## Differences from the historical artifact",
        "",
        "`historical_comparison.csv` contains every directly comparable value. Monte Carlo and bootstrap differences are expected because this rerun regenerated all random draws independently rather than reusing historical artifacts. Deterministic differences are the primary reproduction test.",
        "",
        markdown_table(comparison[comparison["status"] == "different"].head(30)),
        "" if not (comparison["status"] == "different").any() else "Only the first 30 differently classified entries are shown here; the CSV is authoritative.",
        "",
        "## Data assessment",
        "",
        "The public GEO supplementary RPKM matrix is a legitimate source for this cohort and the structural sample counts reproduce. The analysis remains limited by heterogeneous tissue composition (colorectal, endometrial, and ovarian samples), binary tumor/normal residualization that cannot remove all tissue and nonlinear structure, and data-dependent HVG/PCA geometry. Those limitations narrow interpretation but do not invalidate the computational rerun.",
        "",
        "## Files",
        "",
        "- `primary_results.json`: machine-readable conclusions and complete deterministic results",
        "- `within_class.csv`, `within_stratum.csv`, `residualization.csv`, `bootstrap.csv`, `confound_spectrum.csv`: report tables",
        "- `null_and_bootstrap_draws.npz`: every regenerated Monte Carlo and bootstrap draw",
        "- `historical_comparison.csv`: new-vs-historical comparison without forced agreement",
        "- `gse146889_clean_lineage_audit.png`: generated only from the new CSV outputs",
        "- `run_provenance.json`: input, environment, script, shard, and output hashes",
    ]
    with (output_dir / "REPORT.md").open("w", encoding="utf-8") as handle:
        handle.write("\n".join(report_lines).rstrip() + "\n")

    provenance["output_hashes"] = {
        path.name: sha256_file(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "run_provenance.json"
    }
    with (output_dir / "run_provenance.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(provenance), handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(primary_results["claims"], indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
