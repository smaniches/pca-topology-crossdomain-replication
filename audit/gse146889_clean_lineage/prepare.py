#!/usr/bin/env python3
"""Prepare a clean-lineage GSE146889 confound-attribution rerun.

This implementation is independent of the historical transcript-reconstructed
GSE146889 driver and its saved null artifacts. It downloads the public GEO
matrix, fixes the cohort and preprocessing deterministically, computes all
non-Monte-Carlo quantities, and writes a compact prepared artifact for sharded
null/bootstrap workers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
import sklearn
from ripser import ripser
from scipy.stats import fisher_exact, pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146889/"
    "suppl/GSE146889_GeneCount.tsv.gz"
)
EXPECTED_N = 176
EXPECTED_TUMOR = 91
EXPECTED_NORMAL = 85
N_HVG = 2000
N_PCS = 50
PCA_SEED = 42


@dataclass(frozen=True)
class PreparedSubset:
    name: str
    indices: np.ndarray
    standardized: np.ndarray
    pca: np.ndarray
    selected_gene_ids: np.ndarray
    explained_variance_ratio_sum: float
    observed_max_h1: float
    hvg_cutoff_variance: float
    hvg_cutoff_tie_count: int
    hvg_cutoff_selected_tie_count: int
    hvg_boundary_tie_ambiguous: bool


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


def max_h1(point_cloud: np.ndarray) -> float:
    result = ripser(point_cloud, maxdim=1)
    h1 = np.asarray(result["dgms"][1], dtype=float)
    if h1.size == 0:
        return 0.0
    finite = h1[np.isfinite(h1[:, 1])]
    if finite.size == 0:
        return 0.0
    return float(np.max(finite[:, 1] - finite[:, 0]))


def most_persistent_cocycle_support(
    point_cloud: np.ndarray, labels: np.ndarray
) -> dict[str, Any]:
    result = ripser(point_cloud, maxdim=1, do_cocycles=True)
    h1 = np.asarray(result["dgms"][1], dtype=float)
    if h1.size == 0:
        return {
            "persistence": 0.0,
            "support_indices": [],
            "support_n": 0,
            "support_tumor": 0,
            "support_normal": 0,
            "fisher_p_tumor_enrichment": 1.0,
        }
    finite_mask = np.isfinite(h1[:, 1])
    if not np.any(finite_mask):
        return {
            "persistence": 0.0,
            "support_indices": [],
            "support_n": 0,
            "support_tumor": 0,
            "support_normal": 0,
            "fisher_p_tumor_enrichment": 1.0,
        }
    persistence = h1[:, 1] - h1[:, 0]
    persistence[~finite_mask] = -np.inf
    bar_index = int(np.argmax(persistence))
    cocycles_h1 = result.get("cocycles", [[], []])[1]
    if bar_index >= len(cocycles_h1):
        support = np.array([], dtype=int)
    else:
        cocycle = np.asarray(cocycles_h1[bar_index])
        if cocycle.size == 0:
            support = np.array([], dtype=int)
        else:
            support = np.unique(cocycle[:, :2].astype(int).ravel())
    tumor = labels.astype(bool)
    support_mask = np.zeros(len(labels), dtype=bool)
    support_mask[support] = True
    a = int(np.sum(support_mask & tumor))
    b = int(np.sum(support_mask & ~tumor))
    c = int(np.sum(~support_mask & tumor))
    d = int(np.sum(~support_mask & ~tumor))
    p_value = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
    return {
        "bar_index": bar_index,
        "birth": float(h1[bar_index, 0]),
        "death": float(h1[bar_index, 1]),
        "persistence": float(persistence[bar_index]),
        "support_indices": support.tolist(),
        "support_n": int(len(support)),
        "support_tumor": a,
        "support_normal": b,
        "fisher_p_tumor_enrichment": p_value,
        "interpretation_boundary": (
            "Support vertices of the representative ripser H1 cocycle. This is not "
            "a canonical geometric cycle and is used only as a descriptive diagnostic."
        ),
    }


def cv_auc(X: np.ndarray, labels: np.ndarray, *, c_value: float) -> dict[str, Any]:
    y = labels.astype(int)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    model = LogisticRegression(
        C=c_value,
        max_iter=5000,
        solver="lbfgs",
        random_state=42,
    )
    scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
    return {
        "C": c_value,
        "fold_scores": scores.tolist(),
        "mean": float(scores.mean()),
        "sd_population": float(scores.std(ddof=0)),
    }


def auc_permutation_test(
    X: np.ndarray,
    labels: np.ndarray,
    *,
    c_value: float,
    n_perm: int,
    seed: int,
) -> dict[str, Any]:
    observed = cv_auc(X, labels, c_value=c_value)
    rng = np.random.default_rng(seed)
    null_values = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        permuted = rng.permutation(labels)
        null_values[i] = cv_auc(X, permuted, c_value=c_value)["mean"]
    obs_dev = abs(observed["mean"] - 0.5)
    null_dev = np.abs(null_values - 0.5)
    p_value = (1 + int(np.sum(null_dev >= obs_dev))) / (n_perm + 1)
    return {
        "observed": observed,
        "null_mean": float(null_values.mean()),
        "null_sd_population": float(null_values.std(ddof=0)),
        "p_two_sided_distance_from_0_5": float(p_value),
        "n_permutations": n_perm,
        "seed": seed,
    }


def direct_projection_auc(projection: np.ndarray, labels: np.ndarray) -> float:
    auc = float(roc_auc_score(labels.astype(int), projection))
    return max(auc, 1.0 - auc)


def preprocess_subset(
    X_log_all: np.ndarray,
    gene_ids: np.ndarray,
    indices: np.ndarray,
    name: str,
) -> PreparedSubset:
    X = np.asarray(X_log_all[indices], dtype=np.float64)
    variance = X.var(axis=0, ddof=0)
    nonzero = variance > 0
    if int(np.sum(nonzero)) < N_HVG:
        raise RuntimeError(
            f"{name}: only {int(np.sum(nonzero))} nonzero-variance genes; need {N_HVG}"
        )
    X_nz = X[:, nonzero]
    ids_nz = gene_ids[nonzero]
    var_nz = variance[nonzero]

    order = np.argsort(var_nz, kind="mergesort")
    selected_positions = order[-N_HVG:]
    cutoff = float(var_nz[selected_positions].min())
    tied = np.isclose(var_nz, cutoff, rtol=0.0, atol=0.0)
    selected_tied = np.isclose(
        var_nz[selected_positions], cutoff, rtol=0.0, atol=0.0
    )
    tie_count = int(np.sum(tied))
    selected_tie_count = int(np.sum(selected_tied))
    ambiguous = tie_count > selected_tie_count

    X_hvg = X_nz[:, selected_positions]
    selected_ids = ids_nz[selected_positions]
    standardized = StandardScaler().fit_transform(X_hvg)
    if not np.isfinite(standardized).all():
        raise RuntimeError(f"{name}: non-finite standardized matrix")

    n_components = min(N_PCS, standardized.shape[0] - 1, standardized.shape[1])
    pca_model = PCA(n_components=n_components, random_state=PCA_SEED)
    pca = pca_model.fit_transform(standardized)
    observed = max_h1(pca)
    return PreparedSubset(
        name=name,
        indices=np.asarray(indices, dtype=int),
        standardized=standardized,
        pca=pca,
        selected_gene_ids=np.asarray(selected_ids, dtype=str),
        explained_variance_ratio_sum=float(
            np.sum(pca_model.explained_variance_ratio_)
        ),
        observed_max_h1=observed,
        hvg_cutoff_variance=cutoff,
        hvg_cutoff_tie_count=tie_count,
        hvg_cutoff_selected_tie_count=selected_tie_count,
        hvg_boundary_tie_ambiguous=ambiguous,
    )


def residualize_class_means(
    X_standardized: np.ndarray, labels: np.ndarray
) -> tuple[np.ndarray, float]:
    tumor = labels.astype(bool)
    grand = X_standardized.mean(axis=0)
    tumor_mean = X_standardized[tumor].mean(axis=0)
    normal_mean = X_standardized[~tumor].mean(axis=0)
    residualized = X_standardized.copy()
    residualized[tumor] = X_standardized[tumor] - tumor_mean + grand
    residualized[~tumor] = X_standardized[~tumor] - normal_mean + grand
    remaining = float(
        np.max(
            np.abs(
                residualized[tumor].mean(axis=0)
                - residualized[~tumor].mean(axis=0)
            )
        )
    )
    return residualized, remaining


def subset_summary(subset: PreparedSubset, labels: np.ndarray) -> dict[str, Any]:
    local_labels = labels[subset.indices]
    return {
        "name": subset.name,
        "n": int(len(subset.indices)),
        "n_tumor": int(np.sum(local_labels)),
        "n_normal": int(np.sum(~local_labels)),
        "observed_max_h1": subset.observed_max_h1,
        "pca_explained_variance_ratio_sum": subset.explained_variance_ratio_sum,
        "selected_gene_count": int(len(subset.selected_gene_ids)),
        "selected_gene_ids_sha256": hashlib.sha256(
            "\n".join(subset.selected_gene_ids.tolist()).encode("utf-8")
        ).hexdigest(),
        "hvg_cutoff_variance": subset.hvg_cutoff_variance,
        "hvg_cutoff_tie_count": subset.hvg_cutoff_tie_count,
        "hvg_cutoff_selected_tie_count": subset.hvg_cutoff_selected_tie_count,
        "hvg_boundary_tie_ambiguous": subset.hvg_boundary_tie_ambiguous,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--protocol-path", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    data_dir = Path(args.data_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    protocol_path = Path(args.protocol_path)

    data_path = data_dir / "GSE146889_GeneCount.tsv.gz"
    if not data_path.exists():
        print(f"Downloading {GEO_URL}", flush=True)
        urllib.request.urlretrieve(GEO_URL, data_path)
    data_sha = sha256_file(data_path)
    data_size = data_path.stat().st_size
    print(f"Input bytes={data_size} sha256={data_sha}", flush=True)

    frame = pd.read_csv(data_path, sep="\t", low_memory=False)
    rpkm_columns = [column for column in frame.columns if column.endswith("_rpkm")]
    if not rpkm_columns:
        raise RuntimeError("No *_rpkm columns found")
    numeric = frame[rpkm_columns].apply(pd.to_numeric, errors="coerce")
    invalid_rows = numeric.isna().any(axis=1)
    valid_numeric = numeric.loc[~invalid_rows]
    if valid_numeric.empty:
        raise RuntimeError("No valid RPKM rows remain")

    gene_id_column = "GeneId" if "GeneId" in frame.columns else frame.columns[0]
    gene_ids = frame.loc[~invalid_rows, gene_id_column].astype(str).to_numpy()
    rpkm = valid_numeric.to_numpy(dtype=np.float64).T
    if np.any(rpkm < 0):
        raise RuntimeError("Negative RPKM values encountered")

    sample_names = np.asarray([name[: -len("_rpkm")] for name in rpkm_columns])
    tumor_marker = np.char.find(np.char.lower(sample_names), "_tumor_") >= 0
    normal_marker = np.char.find(np.char.lower(sample_names), "_normal_") >= 0
    ambiguous = tumor_marker == normal_marker
    if np.any(ambiguous):
        bad = sample_names[ambiguous].tolist()
        raise RuntimeError(f"Ambiguous tumor/normal labels: {bad[:10]}")
    labels = tumor_marker.astype(bool)
    if len(labels) != EXPECTED_N:
        raise RuntimeError(f"Expected {EXPECTED_N} samples, found {len(labels)}")
    if int(np.sum(labels)) != EXPECTED_TUMOR:
        raise RuntimeError(
            f"Expected {EXPECTED_TUMOR} tumors, found {int(np.sum(labels))}"
        )
    if int(np.sum(~labels)) != EXPECTED_NORMAL:
        raise RuntimeError(
            f"Expected {EXPECTED_NORMAL} normals, found {int(np.sum(~labels))}"
        )

    X_log_all = np.log1p(rpkm)
    all_indices = np.arange(len(labels), dtype=int)
    mixed = preprocess_subset(X_log_all, gene_ids, all_indices, "mixed")
    tumor = preprocess_subset(
        X_log_all, gene_ids, np.flatnonzero(labels), "tumor"
    )
    normal = preprocess_subset(
        X_log_all, gene_ids, np.flatnonzero(~labels), "normal"
    )

    mean_shift = mixed.standardized[labels].mean(axis=0) - mixed.standardized[
        ~labels
    ].mean(axis=0)
    norm = float(np.linalg.norm(mean_shift))
    if norm == 0.0:
        raise RuntimeError("Zero class-mean-shift vector")
    direction = mean_shift / norm
    projection = mixed.standardized @ direction
    order = np.argsort(projection, kind="mergesort")
    quartile_indices = [np.sort(chunk) for chunk in np.array_split(order, 4)]
    if [len(chunk) for chunk in quartile_indices] != [44, 44, 44, 44]:
        raise RuntimeError("Quartile sizes are not 44 each")
    quartiles = [
        preprocess_subset(X_log_all, gene_ids, indices, f"q{i + 1}")
        for i, indices in enumerate(quartile_indices)
    ]

    residualized_std, class_mean_diff = residualize_class_means(
        mixed.standardized, labels
    )
    residualized_pca_model = PCA(n_components=N_PCS, random_state=PCA_SEED)
    residualized_pca = residualized_pca_model.fit_transform(residualized_std)
    residualized_observed = max_h1(residualized_pca)

    intact_cocycle = most_persistent_cocycle_support(mixed.pca, labels)
    residualized_cocycle = most_persistent_cocycle_support(
        residualized_pca, labels
    )

    auc_intact_c001 = cv_auc(mixed.pca, labels, c_value=0.01)
    auc_resid_c001 = cv_auc(residualized_pca, labels, c_value=0.01)
    auc_intact_c1 = cv_auc(mixed.pca, labels, c_value=1.0)
    auc_resid_c1 = cv_auc(residualized_pca, labels, c_value=1.0)
    auc_perm = auc_permutation_test(
        residualized_pca,
        labels,
        c_value=0.01,
        n_perm=200,
        seed=14689300,
    )

    pc_correlations: list[dict[str, Any]] = []
    for index in range(mixed.pca.shape[1]):
        correlation, p_value = pearsonr(mixed.pca[:, index], projection)
        pc_correlations.append(
            {
                "pc": index + 1,
                "pearson_r": float(correlation),
                "pearson_p": float(p_value),
            }
        )
    variance_total = float(
        np.sum(np.var(mixed.standardized, axis=0, ddof=0))
    )
    variance_direction = float(np.var(projection, ddof=0))

    tie_sensitivity: dict[str, Any] = {"performed": False}
    if mixed.hvg_boundary_tie_ambiguous:
        full_variance = X_log_all.var(axis=0, ddof=0)
        nonzero = full_variance > 0
        X_nz = X_log_all[:, nonzero]
        ids_nz = gene_ids[nonzero]
        var_nz = full_variance[nonzero]
        legacy_positions = np.argsort(var_nz)[-N_HVG:]
        legacy_std = StandardScaler().fit_transform(X_nz[:, legacy_positions])
        legacy_pca = PCA(n_components=N_PCS, random_state=PCA_SEED).fit_transform(
            legacy_std
        )
        tie_sensitivity = {
            "performed": True,
            "legacy_selected_gene_ids_sha256": hashlib.sha256(
                "\n".join(ids_nz[legacy_positions].astype(str).tolist()).encode(
                    "utf-8"
                )
            ).hexdigest(),
            "legacy_observed_max_h1": max_h1(legacy_pca),
            "stable_observed_max_h1": mixed.observed_max_h1,
        }

    real_results = {
        "data": {
            "url": GEO_URL,
            "path": str(data_path),
            "sha256": data_sha,
            "size_bytes": data_size,
            "table_rows_total": int(len(frame)),
            "invalid_or_missing_rpkm_rows_removed": int(np.sum(invalid_rows)),
            "valid_gene_rows": int(valid_numeric.shape[0]),
            "rpkm_column_count": len(rpkm_columns),
            "sample_count": len(labels),
            "tumor_count": int(np.sum(labels)),
            "normal_count": int(np.sum(~labels)),
            "label_rule": "sample name contains exactly one of _tumor_ or _normal_",
        },
        "preprocessing": {
            "transform": "log1p(RPKM)",
            "variance_ddof": 0,
            "hvg_count": N_HVG,
            "hvg_sort": "numpy.argsort(kind='mergesort'), highest variances",
            "standardization": "sklearn.preprocessing.StandardScaler per analysis subset",
            "pca_components": N_PCS,
            "pca_random_state": PCA_SEED,
            "persistent_homology": "ripser maxdim=1; maximum finite H1 persistence",
        },
        "subsets": {
            subset.name: subset_summary(subset, labels)
            for subset in [mixed, tumor, normal, *quartiles]
        },
        "quartiles": [
            {
                "name": subset.name,
                "global_indices": subset.indices.tolist(),
                "projection_min": float(projection[subset.indices].min()),
                "projection_max": float(projection[subset.indices].max()),
                "n_tumor": int(np.sum(labels[subset.indices])),
                "n_normal": int(np.sum(~labels[subset.indices])),
            }
            for subset in quartiles
        ],
        "residualization": {
            "observed_max_h1": residualized_observed,
            "pca_explained_variance_ratio_sum": float(
                np.sum(residualized_pca_model.explained_variance_ratio_)
            ),
            "max_abs_class_mean_difference_after": class_mean_diff,
            "auc_C_0_01_intact": auc_intact_c001,
            "auc_C_0_01_residualized": auc_resid_c001,
            "auc_C_1_intact": auc_intact_c1,
            "auc_C_1_residualized": auc_resid_c1,
            "auc_C_0_01_residualized_permutation_test": auc_perm,
            "intact_cocycle_support": intact_cocycle,
            "residualized_cocycle_support": residualized_cocycle,
        },
        "confound_spectrum": {
            "direction_norm_before_normalization": norm,
            "projection_direct_auc_orientation_invariant": direct_projection_auc(
                projection, labels
            ),
            "variance_along_direction": variance_direction,
            "total_standardized_feature_variance": variance_total,
            "variance_fraction_along_direction": variance_direction / variance_total,
            "pc_correlations": pc_correlations,
        },
        "hvg_tie_sensitivity": tie_sensitivity,
    }

    np.savez_compressed(
        output_dir / "prepared.npz",
        labels=labels.astype(np.uint8),
        sample_names=sample_names.astype(str),
        mixed_std=mixed.standardized,
        mixed_pca=mixed.pca,
        tumor_std=tumor.standardized,
        tumor_pca=tumor.pca,
        normal_std=normal.standardized,
        normal_pca=normal.pca,
        q1_std=quartiles[0].standardized,
        q1_pca=quartiles[0].pca,
        q2_std=quartiles[1].standardized,
        q2_pca=quartiles[1].pca,
        q3_std=quartiles[2].standardized,
        q3_pca=quartiles[2].pca,
        q4_std=quartiles[3].standardized,
        q4_pca=quartiles[3].pca,
        q1_indices=quartiles[0].indices,
        q2_indices=quartiles[1].indices,
        q3_indices=quartiles[2].indices,
        q4_indices=quartiles[3].indices,
        residualized_std=residualized_std,
        residualized_pca=residualized_pca,
        confound_projection=projection,
    )
    with (output_dir / "real_results.json").open("w", encoding="utf-8") as handle:
        json.dump(jsonable(real_results), handle, indent=2, sort_keys=True)
        handle.write("\n")

    provenance = {
        "protocol_path": str(protocol_path),
        "protocol_sha256": sha256_file(protocol_path),
        "execution_git_sha": os.environ.get("GITHUB_SHA"),
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "workflow_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "input_sha256": data_sha,
        "input_size_bytes": data_size,
        "prepare_script_sha256": sha256_file(Path(__file__)),
    }
    with (output_dir / "provenance_base.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(jsonable(provenance), handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(jsonable(real_results["subsets"]), indent=2), flush=True)
    print(
        f"Residualized observed max-H1={residualized_observed:.6f}; "
        f"class-mean residual={class_mean_diff:.3e}",
        flush=True,
    )


if __name__ == "__main__":
    main()
