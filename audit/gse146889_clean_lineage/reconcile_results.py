#!/usr/bin/env python3
"""Reconcile comparison semantics for the GSE146889 clean-lineage audit.

This script does not rerun any experiment. It operates only on the committed
clean-lineage outputs. It corrects two derived comparison mismatches:

1. The historical AUC was computed with LogisticRegression(C=1.0), while the
   first aggregate report compared it with the newly added C=0.01 sensitivity.
2. The historical loop-enrichment p-value was two-sided Fisher exact, while the
   first aggregate report compared it with the one-sided enrichment p-value.

Both definitions were already computed from the same clean rerun. The raw data,
PCA/PH statistics, null draws, bootstrap draws, and cocycle supports are unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

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


def markdown_table(frame: pd.DataFrame, precision: int = 4) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(map(str, columns)) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        cells: list[str] = []
        for value in row.tolist():
            if isinstance(value, (bool, np.bool_)):
                cells.append("yes" if bool(value) else "no")
            elif isinstance(value, (float, np.floating)):
                if np.isnan(value):
                    cells.append("NaN")
                elif value != 0 and abs(value) < 1e-3:
                    cells.append(f"{value:.{precision}e}")
                else:
                    cells.append(f"{value:.{precision}f}")
            else:
                cells.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
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
        h = float(historical)
        n = float(new)
        absolute = abs(n - h)
        relative = absolute / abs(h) if h != 0 else np.nan
        if absolute <= 1e-9:
            status = "exact"
        elif metric_type == "deterministic" and absolute <= 1e-3:
            status = "numerically_reproduced"
        elif metric_type == "deterministic" and relative <= 0.01:
            status = "close_deterministic"
        elif metric_type in {"monte_carlo", "bootstrap"} and relative <= 0.10:
            status = "within_10pct"
        else:
            status = "different"
    except (TypeError, ValueError):
        h, n = historical, new
        absolute = relative = np.nan
        status = "same" if str(historical) == str(new) else "different"
    rows.append(
        {
            "section": section,
            "key": key,
            "metric": metric,
            "metric_type": metric_type,
            "historical": h,
            "new": n,
            "absolute_difference": absolute,
            "relative_difference": relative,
            "status": status,
        }
    )


def make_comparison(repo_root: Path, results_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    specifications = {
        "within_class": (
            "condition",
            [
                ("observed_max_H1", "deterministic"),
                ("z_gaussian", "monte_carlo"),
                ("p_gaussian", "monte_carlo"),
                ("z_permutation", "monte_carlo"),
                ("p_permutation", "monte_carlo"),
            ],
            results_dir / "within_class.csv",
        ),
        "within_stratum": (
            "quartile",
            [
                ("observed_max_H1", "deterministic"),
                ("z_gauss", "monte_carlo"),
                ("z_perm", "monte_carlo"),
            ],
            results_dir / "within_stratum.csv",
        ),
        "residualization": (
            "space",
            [
                ("observed_max_H1", "deterministic"),
                ("z_gaussian", "monte_carlo"),
                ("z_permutation", "monte_carlo"),
                ("cv_auc_tumor_normal", "deterministic"),
                ("top_loop_n_participants", "diagnostic"),
                ("top_loop_n_tumor", "diagnostic"),
                ("top_loop_n_normal", "diagnostic"),
                ("fisher_p_tumor_enrichment", "diagnostic"),
            ],
            results_dir / "residualization.csv",
        ),
        "bootstrap": (
            "subset",
            [
                ("point_delta", "bootstrap"),
                ("delta_ci95_lo", "bootstrap"),
                ("delta_ci95_hi", "bootstrap"),
                ("point_z", "bootstrap"),
                ("z_ci95_lo", "bootstrap"),
                ("z_ci95_hi", "bootstrap"),
                ("excludes_zero", "bootstrap"),
            ],
            results_dir / "bootstrap.csv",
        ),
    }
    for section, (key_column, metrics, new_path) in specifications.items():
        historical = pd.read_csv(repo_root / HISTORICAL_FILES[section])
        new = pd.read_csv(new_path)
        hist_map = {str(row[key_column]): row for _, row in historical.iterrows()}
        new_map = {str(row[key_column]): row for _, row in new.iterrows()}
        for key in sorted(set(hist_map) & set(new_map)):
            for metric, metric_type in metrics:
                compare_value(
                    rows,
                    section=section,
                    key=key,
                    metric=metric,
                    historical=hist_map[key][metric],
                    new=new_map[key][metric],
                    metric_type=metric_type,
                )
    return pd.DataFrame(rows)


def correct_residualization(primary: dict[str, Any], results_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(results_dir / "residualization.csv")
    real = primary["real_results"]["residualization"]
    mapping = {
        "Full HVG->PCA50 (class-mean intact)": (
            real["auc_C_1_intact"],
            real["auc_C_0_01_intact"],
            real["intact_cocycle_support"],
        ),
        "Class-mean-residualized HVG->PCA50": (
            real["auc_C_1_residualized"],
            real["auc_C_0_01_residualized"],
            real["residualized_cocycle_support"],
        ),
    }
    total_tumor = int(primary["data"]["tumor_count"])
    total_normal = int(primary["data"]["normal_count"])
    for index, row in frame.iterrows():
        c1, c001, support = mapping[str(row["space"])]
        a = int(support["support_tumor"])
        b = int(support["support_normal"])
        c = total_tumor - a
        d = total_normal - b
        p_two = float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
        p_greater = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        frame.loc[index, "cv_auc_tumor_normal"] = float(c1["mean"])
        frame.loc[index, "cv_auc_sd"] = float(c1["sd_population"])
        frame.loc[index, "cv_auc_C_0_01_sensitivity"] = float(c001["mean"])
        frame.loc[index, "cv_auc_C_0_01_sd"] = float(c001["sd_population"])
        frame.loc[index, "fisher_p_tumor_enrichment"] = p_two
        frame.loc[index, "fisher_p_tumor_enrichment_one_sided_sensitivity"] = p_greater
        frame.loc[index, "auc_historical_comparator"] = "LogisticRegression(C=1.0), 5-fold StratifiedKFold"
        frame.loc[index, "fisher_historical_comparator"] = "Fisher exact, two-sided"
    frame.to_csv(results_dir / "residualization.csv", index=False)
    return frame


def regenerate_figure(results_dir: Path) -> None:
    """Write a compact, dependency-light SVG from the corrected CSV outputs."""
    from html import escape

    within_class = pd.read_csv(results_dir / "within_class.csv")
    within_stratum = pd.read_csv(results_dir / "within_stratum.csv")
    residualization = pd.read_csv(results_dir / "residualization.csv")
    bootstrap = pd.read_csv(results_dir / "bootstrap.csv")

    width, height = 1200, 820
    panel_w, panel_h = 520, 310
    positions = [(50, 80), (630, 80), (50, 450), (630, 450)]
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.title{font-size:23px;font-weight:700}.panel{font-size:17px;font-weight:700}.label{font-size:12px}.small{font-size:11px}.axis{stroke:#222;stroke-width:1}.mark{fill:#777}.mark2{fill:#bbb}.line{stroke:#222;stroke-width:2;fill:none}.dash{stroke:#555;stroke-width:1;stroke-dasharray:5 4}.ci{stroke:#222;stroke-width:2}</style>',
        '<text x="600" y="35" text-anchor="middle" class="title">GSE146889 clean-lineage confound-attribution rerun</text>',
        '<text x="600" y="58" text-anchor="middle" class="small">Matched historical definitions; regenerated null and bootstrap draws</text>',
    ]

    def panel_frame(x0: int, y0: int, title: str) -> tuple[int, int, int, int]:
        svg.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="none" stroke="#bbb"/>')
        svg.append(f'<text x="{x0+15}" y="{y0+25}" class="panel">{escape(title)}</text>')
        left, top, right, bottom = x0 + 55, y0 + 50, x0 + panel_w - 20, y0 + panel_h - 55
        svg.append(f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/>')
        svg.append(f'