"""Prospective CPTAC-CCRCC conditional-permutation sensitivity experiment.

Not part of the original preregistration. Preserves historical result artifacts.
Run from repository root:
  python audit/cptac_resid_null/experiment.py --mode pilot --draws 24 --out /tmp/cptac-pilot.json

Null: independently permute each protein *within each class* before applying
the identical class-mean residualization and PCA50 pipeline to every draw.
This preserves univariate class-conditional empirical distributions, not
cross-protein correlations or matched-patient structure. Interpret p-values
only under conditional sample-exchangeability assumptions.
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code" / "confound_attribution_audit"))
from confound_audit_common import max_h1_persistence, residualize_class_mean  # noqa: E402
from confound_audit_cptac_ccrcc import preprocess  # noqa: E402

RAW = ROOT / "code/replication_CPTAC_CCRCC/data/cptac_ccrcc_log2ratio_raw.csv"
LABELS = ROOT / "code/replication_CPTAC_CCRCC/data/cptac_ccrcc_labels.csv"
HASHES = {
    RAW: "d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3",
    LABELS: "5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e",
}


def verify_inputs():
    for path, expected in HASHES.items():
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Input SHA256 mismatch: {path}; {actual} != {expected}")


def prepare():
    verify_inputs()
    raw = pd.read_csv(RAW, index_col=0)
    labels = pd.read_csv(LABELS, index_col=0)
    if not raw.index.is_unique or not labels.index.is_unique:
        raise ValueError("Non-unique sample identifiers")
    if set(raw.index) != set(labels.index):
        raise ValueError("Matrix/label IDs do not match exactly")
    if not set(labels.iloc[:, 0].unique()) <= {"Primary Tumor", "Solid Tissue Normal"}:
        raise ValueError("Unexpected class labels")
    X, mask, normals = preprocess(raw, labels)
    if X.shape != (194, 2000) or (int(mask.sum()), int(normals.sum())) != (110, 84):
        raise ValueError(f"Unexpected cohort shape or class balance: {X.shape}")
    if not np.isfinite(X).all():
        raise ValueError("Nonfinite standardized feature")
    return X, mask


def shuffle_columns_within_classes(X, mask, rng):
    """Conditional univariate-margin-preserving surrogate, with fixed labels."""
    X = np.asarray(X, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if X.ndim != 2 or mask.shape != (X.shape[0],) or np.unique(mask).size != 2:
        raise ValueError("Expected a finite matrix and a binary nonempty class mask")
    if not np.isfinite(X).all():
        raise ValueError("Nonfinite feature")
    out = np.empty_like(X)
    for group in (False, True):
        ix = np.flatnonzero(mask == group)
        for j in range(X.shape[1]):
            out[ix, j] = rng.permutation(X[ix, j])
    return out


def residualized_pca(X, mask, n_pcs=50):
    Xr, max_class_mean_diff = residualize_class_mean(X, mask)
    if max_class_mean_diff > 1e-9:
        raise ValueError(f"Class means not removed: {max_class_mean_diff}")
    components = min(n_pcs, X.shape[0] - 1, X.shape[1])
    # Same estimator and seed for observed and every surrogate.
    return PCA(n_components=components, random_state=42).fit_transform(Xr)


def statistic(X, mask):
    # Ripser is applied to *both* observed and null draws, not compared with
    # Gudhi-generated nulls from the historical audit.
    return max_h1_persistence(residualized_pca(X, mask), backend="ripser")


def experiment(X, mask, draws, seed, max_seconds=1500):
    if draws < 2 or max_seconds <= 0:
        raise ValueError("Need >=2 null draws and positive time limit")
    t0 = time.monotonic()
    observed = statistic(X, mask)
    # Backend / preprocessing guard against a silent drift from historical real value.
    if abs(observed - 4.409938425199396) > 0.08:
        raise ValueError(f"Unexpected observed H1={observed:.6f}, historical Gudhi=4.409938")
    rng = np.random.default_rng(seed)
    null = []
    for i in range(draws):
        if time.monotonic() - t0 > max_seconds:
            raise TimeoutError(f"Time cap reached after {i}/{draws} null draws")
        sur = shuffle_columns_within_classes(X, mask, rng)
        v = statistic(sur, mask)
        if not np.isfinite(v):
            raise ValueError(f"Nonfinite null draw {i}")
        null.append(float(v))
        if (i + 1) % 8 == 0:
            print(f"completed {i+1}/{draws} draws", flush=True)
    null = np.asarray(null)
    sd = float(null.std(ddof=1))
    return {
        "observed_max_h1": float(observed),
        "null_max_h1": null.tolist(),
        "null_mean": float(null.mean()),
        "null_sd": sd,
        "z_descriptive": float((observed - null.mean()) / sd) if sd > 0 else None,
        "p_one_sided_plus_one": float((1 + (null >= observed).sum()) / (len(null) + 1)),
        "num_null_ge_observed": int((null >= observed).sum()),
        "draws": draws,
        "seed": seed,
        "seconds": round(time.monotonic() - t0, 2),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("pilot", "confirm"), required=True)
    parser.add_argument("--draws", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    draws = args.draws if args.draws is not None else (24 if args.mode == "pilot" else 499)
    if args.mode == "pilot" and not (2 <= draws <= 32):
        parser.error("pilot permits 2-32 draws")
    if args.mode == "confirm" and draws != 499:
        parser.error("confirmation requires exactly 499 draws")
    X, mask = prepare()
    payload = experiment(X, mask, draws=draws, seed=args.seed,
                         max_seconds=1500 if args.mode == "pilot" else 6900)
    payload.update({
        "mode": args.mode,
        "status": "EXPLORATORY" if args.mode == "pilot" else "PROSPECTIVE_SENSITIVITY",
        "null_interpretation": "Within-label, per-column shuffle; presumes conditional sample exchangeability; does not preserve cross-feature covariance or patient pairing",
        "data_sha256": {p.relative_to(ROOT).as_posix(): h for p, h in HASHES.items()},
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "sklearn": sklearn.__version__,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out}: p={payload['p_one_sided_plus_one']:.6f}; "
          f"z={payload['z_descriptive']}; elapsed={payload['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
