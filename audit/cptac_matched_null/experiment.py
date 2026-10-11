"""CPTAC-TUMOR-MATCHED-NULL-002: tumor-only featurewise permutation control.

Post-preregistration reference model, frozen in PROTOCOL.md. Only the null
construction differs from the archived Haar covariance-preserving analysis:
same prepared 110 x 2000 tumor matrix, PCA(50), Euclidean ripser max finite H1.

This software does NOT establish biological causation, batch independence,
patient exchangeability, or a null-model-invariant topological mechanism.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.linalg import eigh
from scipy.spatial.distance import pdist
from scipy.stats import beta
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
HAAR_JSON = ROOT / "results/cptac_tumor_covariance_null_20261009/499_draws.json"
HAAR_JSON_SHA256 = "e45385676be011c6f2df7f6dbbe30202d5a28418bc6e490e8af176e77b89fcfa"
EXPECTED_OBSERVED = 5.354297637939453
K = 50
DRAW_COUNT = 499
PILOT_SEED = 20261012
PRIMARY_SEED = 20261011
ID = "CPTAC-TUMOR-MATCHED-NULL-002"


def haar_source():
    """Import exact observed-cloud preparation and ripser statistic from prior audit."""
    source = ROOT / "audit/cptac_covariance_null/experiment.py"
    spec = importlib.util.spec_from_file_location("cptac_frozen_haar_impl", source)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import frozen covariance-null source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def archived_reference():
    """Fail closed on unexpected modification of all archived Haar draws."""
    digest = sha256_bytes(HAAR_JSON.read_bytes())
    if digest != HAAR_JSON_SHA256:
        raise ValueError("Archived covariance-preserving null was changed: " + digest)
    d = json.loads(HAAR_JSON.read_text(encoding="utf-8"))
    v = np.asarray(d["null_max_h1"], dtype=np.float64)
    if (len(v) != DRAW_COUNT or d["draws"] != DRAW_COUNT or
            d["null_exceedances"] != 45 or
            not np.isclose(d["observed_max_h1"], EXPECTED_OBSERVED, atol=1e-9) or
            np.count_nonzero(v >= EXPECTED_OBSERVED) != 45 or
            not np.isclose(d["p_one_sided_plus_one"], 0.092, atol=1e-14)):
        raise ValueError("Archived covariance-null statistics no longer match evidence")
    return d


def permute_features_independently(x, rng):
    """Each protein's values are permuted across the same tumor rows."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] < 4 or not np.isfinite(x).all():
        raise ValueError("Expected finite 2D sample-by-protein matrix")
    n, p = x.shape
    out = np.empty_like(x)
    for j in range(p):
        out[:, j] = x[rng.permutation(n), j]
    return out


def exact_pca_scores(x, components=K):
    """Efficient full-rank PCA scores from eigenvectors of centered row Gram.

    Same pairwise Euclidean distances as sklearn PCA with full SVD, after
    independent parity checks. Avoids repeated wide 110-by-2000 SVDs.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2 or not np.isfinite(x).all():
        raise ValueError("Expected finite matrix")
    n, p = x.shape
    if not 1 <= components <= min(n - 1, p):
        raise ValueError("Invalid PCA dimensionality")
    centered = x - x.mean(axis=0)
    gram = centered @ centered.T
    eigenvalues, eigenvectors = eigh(
        gram, subset_by_index=(n-components, n-1), check_finite=True
    )
    eig = eigenvalues[::-1]
    if np.any(eig <= 0) or not np.isfinite(eig).all():
        raise ValueError("Nonpositive retained PCA eigenvalue")
    scores = eigenvectors[:, ::-1] * np.sqrt(eig)[None, :]
    if not np.isfinite(scores).all():
        raise ValueError("Nonfinite PCA scores")
    return scores


def pca_distance_parity(x, scores, max_error=1e-8):
    """Independently verify PCA point-cloud metric using full sklearn SVD."""
    direct = PCA(n_components=scores.shape[1], svd_solver="full").fit_transform(x)
    discrepancy = float(np.max(np.abs(pdist(scores) - pdist(direct))))
    if discrepancy > max_error or not np.isfinite(discrepancy):
        raise ValueError("PCA Gram / full-SVD metric disagrees: " + str(discrepancy))
    return discrepancy


def clopper_pearson(k, n, alpha=0.05):
    """Interval for Monte Carlo exceedance fraction, not a biological effect."""
    lo = float(beta.ppf(alpha/2, k, n-k+1)) if k else 0.0
    hi = float(beta.ppf(1-alpha/2, k+1, n-k)) if k < n else 1.0
    return [lo, hi]


def run(x, draws, seed, seconds_limit=1500):
    if draws < 2 or seconds_limit <= 0:
        raise ValueError("Need >=2 null draws and a positive execution budget")
    old = archived_reference()
    source = haar_source()
    x = np.asarray(x, dtype=np.float64)
    if x.shape != (110, 2000) or not np.isfinite(x).all():
        raise ValueError("Observed tumor matrix changed shape or finiteness")
    started = time.monotonic()
    observed_scores = source.score_state(x, K)["scores"]
    observed = float(source.h1_max(observed_scores))
    if abs(observed - EXPECTED_OBSERVED) > 2e-5:
        raise ValueError("Different observed tumor H1: " + str(observed))
    observed_distance_err = pca_distance_parity(x, observed_scores)
    rng = np.random.default_rng(seed)
    null = np.empty(draws, dtype=np.float64)
    first_permutation_distance_error = None

    for index in range(draws):
        if time.monotonic() - started > seconds_limit:
            raise TimeoutError(f"Execution time exceeded at {index}/{draws}")
        shuffled = permute_features_independently(x, rng)
        if index == 0:
            # These checks establish that the null preserves all single-protein
            # empirical distributions without silently changing the data space.
            if not np.array_equal(np.sort(shuffled, axis=0), np.sort(x, axis=0)):
                raise ValueError("Independent feature null altered protein marginals")
            if np.array_equal(shuffled, x):
                raise ValueError("Unexpected identity feature permutation")
        pca_scores = exact_pca_scores(shuffled, K)
        if index == 0:
            first_permutation_distance_error = pca_distance_parity(shuffled, pca_scores)
        null[index] = float(source.h1_max(pca_scores))
        if not np.isfinite(null[index]) or null[index] < 0:
            raise ValueError("Invalid null topology")
        if (index + 1) % 25 == 0:
            print(f"completed {index+1}/{draws} matched feature-permutation nulls",flush=True)

    exceedances = int(np.count_nonzero(null >= observed))
    mean = float(np.mean(null))
    sd = float(np.std(null, ddof=1))
    old_vals = np.asarray(old["null_max_h1"], dtype=np.float64)
    haar_exceedances = int(np.count_nonzero(old_vals >= observed))
    if haar_exceedances != 45:
        raise ValueError("Archived Haar exceedance count drift")
    return {
        "experiment_id": ID,
        "status": "POST_PREREGISTRATION_OUTCOME_INFORMED_SENSITIVITY",
        "observed_max_h1": observed,
        "sample_count": x.shape[0],
        "features": x.shape[1],
        "pca_components": K,
        "metric": "Euclidean",
        "topology_statistic": "longest finite ripser H1 persistence",
        "backend": "ripser",
        "new_reference": {
            "type": "independent within-tumor feature-column permutations",
            "null_max_h1": null.tolist(),
            "null_mean": mean,
            "null_sd": sd,
            "null_min": float(null.min()),
            "null_max": float(null.max()),
            "draws": draws,
            "seed": seed,
            "exceedances": exceedances,
            "p_one_sided_plus_one": (1+exceedances)/(draws+1),
            "monte_carlo_exceedance_fraction_95pct_interval":
                clopper_pearson(exceedances, draws),
            "descriptive_z": (observed-mean)/sd if sd>0 else None,
        },
        "archived_covariance_preserving_reference": {
            "json_sha256": HAAR_JSON_SHA256,
            "null_mean": float(old_vals.mean()),
            "null_sd": float(old_vals.std(ddof=1)),
            "draws": len(old_vals),
            "seed": int(old["seed"]),
            "exceedances": haar_exceedances,
            "p_one_sided_plus_one": (haar_exceedances+1)/(len(old_vals)+1),
            "monte_carlo_exceedance_fraction_95pct_interval":
                clopper_pearson(haar_exceedances, len(old_vals)),
        },
        "checks": {
            "observed_metric_full_svd_max_abs_error": observed_distance_err,
            "first_null_metric_full_svd_max_abs_error":
                first_permutation_distance_error,
            "input_sha256": source.SOURCE_HASHES,
        },
        "limitations": [
            "Two null families retain different invariants; contrasting their p-values does not prove covariance causation.",
            "Patient pairing, batch structure, sample exchangeability and feature-selection conditioning remain unverified.",
            "Both nulls use the same tumor cohort but neither is an independent biological replication.",
            "The earlier mixed tumor-normal class-residualized p=0.002 is not compared here.",
        ],
        "elapsed_seconds": round(time.monotonic()-started, 4),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=["pilot","confirm"])
    parser.add_argument("--draws", type=int, default=None)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    count = args.draws if args.draws is not None else (24 if args.mode=="pilot" else DRAW_COUNT)
    if args.mode=="pilot" and not 2 <= count <= 32:
        parser.error("Pilot is bounded to 2..32 null draws")
    if args.mode=="confirm" and count != DRAW_COUNT:
        parser.error("Primary comparison requires exactly 499 null draws")
    source = haar_source()
    x = source.load_tumor()  # verifies cached PDC input hashes/labels
    seed = PILOT_SEED if args.mode=="pilot" else PRIMARY_SEED
    payload = run(x, count, seed, 600 if args.mode=="pilot" else 5400)
    payload.update({
        "mode": args.mode,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "sklearn": sklearn.__version__,
        "source_commit": os.environ.get("GITHUB_SHA", ""),
        "workflow_run": os.environ.get("GITHUB_RUN_ID", ""),
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    current = payload["new_reference"]
    print(f"MATCHED_NULL_RESULT: obs={payload['observed_max_h1']:.9f} "
          f"feature-p={current['p_one_sided_plus_one']:.6f} "
          f"exceedances={current['exceedances']}/{current['draws']}; "
          f"haar-p={payload['archived_covariance_preserving_reference']['p_one_sided_plus_one']:.6f}; "
          f"seconds={payload['elapsed_seconds']}",flush=True)


if __name__ == "__main__":
    main()
