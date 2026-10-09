"""CPTAC-TUMOR-COV-NULL-001: a post-preregistration conditional covariance null.

The exact mean-centered feature scatter matrix is invariant under the Haar
sample-space rotations. Full-SVD PCA scores are generated algebraically from
one source SVD (not estimated independently on a permuted feature table).

This tests a conditional iid matrix-normal row model, not biological causality.
See PROTOCOL.md before interpreting a tail probability.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import sklearn
from scipy.linalg import svd
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
SOURCE_HASHES = {
    "code/replication_CPTAC_CCRCC/data/cptac_ccrcc_log2ratio_raw.csv":
        "d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3",
    "code/replication_CPTAC_CCRCC/data/cptac_ccrcc_labels.csv":
        "5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e",
}
K = 50
H1_BACKEND = "ripser"


def helmert_basis(n):
    """Fixed orthonormal basis of the sample-mean-zero subspace."""
    if n < 4:
        raise ValueError("At least four independent rows are required")
    u = np.zeros((n, n - 1), dtype=np.float64)
    for j in range(n - 1):
        denom = np.sqrt((j + 1) * (j + 2))
        u[:j + 1, j] = 1.0 / denom
        u[j + 1, j] = -(j + 1) / denom
    return u


def haar_orthogonal(d, rng):
    """Haar orthogonal group sample; sign correction removes QR bias."""
    if d < 1:
        raise ValueError("Expected positive orthogonal-group dimension")
    a = rng.standard_normal((d, d))
    q, r = np.linalg.qr(a)
    signs = np.where(np.diag(r) < 0.0, -1.0, 1.0)
    return q * signs[np.newaxis, :]


def score_state(x, components=K):
    """Exact PCA coordinates from the thin SVD in Helmert sample space."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2 or not np.all(np.isfinite(x)):
        raise ValueError("Expected a finite 2-D sample-by-feature matrix")
    n, p = x.shape
    if not (1 <= components <= min(n - 1, p)):
        raise ValueError("PCA component count exceeds sample-centered rank bound")
    u = helmert_basis(n)
    mu = x.mean(axis=0)
    z = u.T @ (x - mu)
    left, singular, _ = svd(z, full_matrices=False, check_finite=True)
    scores = (u @ left[:, :components]) * singular[np.newaxis, :components]
    if not np.all(np.isfinite(scores)):
        raise ValueError("PCA score calculation returned nonfinite numbers")
    return {
        "u": u, "mu": mu, "z": z, "left": left[:, :components],
        "singular": singular, "scores": scores,
        "n": n, "p": p, "components": components,
    }


def rotated_scores(state, rng, check=False):
    """Draw exact covariance-preserving surrogate PCA scores."""
    q = haar_orthogonal(state["n"] - 1, rng)
    rotated = (state["u"] @ (q @ state["left"])) * state["singular"][np.newaxis, :state["components"]]
    if check:
        if np.max(np.abs(q.T @ q - np.eye(len(q)))) > 1e-10:
            raise ValueError("Haar matrix failed orthogonality")
        # Check actual covariance, not just a second-order surrogate scalar.
        # The small feature submatrix avoids a quadratic allocation in 2,000 proteins.
        idx = np.linspace(0, state["p"] - 1, min(17, state["p"]), dtype=int)
        zx = state["z"][:, idx]
        zrot = q @ zx
        if np.max(np.abs(zrot.T @ zrot - zx.T @ zx)) > 1e-7:
            raise ValueError("Surrogate does not preserve covariance")
        if np.max(np.abs((state["u"] @ (q @ zx)).mean(axis=0))) > 1e-10:
            raise ValueError("Surrogate changed feature means")
    return rotated


def h1_max(scores):
    from ripser import ripser
    diagram = np.asarray(ripser(scores, maxdim=1)["dgms"][1], dtype=float)
    if diagram.size == 0:
        return 0.0
    finite = diagram[np.isfinite(diagram[:, 1])]
    if finite.size == 0:
        return 0.0
    out = float(np.max(finite[:, 1] - finite[:, 0]))
    if not np.isfinite(out) or out < 0:
        raise ValueError("Invalid H1 lifetime")
    return out


def verify_inputs():
    for rel, expected in SOURCE_HASHES.items():
        data = ROOT / rel
        actual = hashlib.sha256(data.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Data hash mismatch for {rel}: {actual}")


def load_tumor():
    verify_inputs()
    sys.path.insert(0, str(ROOT / "audit" / "cptac_resid_null"))
    from experiment import prepare  # existing hashed cohort, normalization and label checks
    x, tumor_mask = prepare()
    tumor = x[tumor_mask]
    if tumor.shape != (110, 2000):
        raise ValueError(f"Unexpected tumor subset {tumor.shape}")
    return tumor


def verify_pca_parity(x, state, tol=1e-6):
    """Independent full-SVD reference to detect a wiring or shape error."""
    independent = PCA(n_components=state["components"], svd_solver="full").fit_transform(x)
    err = float(np.max(np.abs(pdist(independent) - pdist(state["scores"]))))
    if not np.isfinite(err) or err > tol:
        raise ValueError(f"Full-SVD PCA pairwise-distance mismatch: {err}")
    return err


def experiment(x, draws, seed, max_seconds):
    if draws < 2:
        raise ValueError("At least two null draws required")
    start = time.monotonic()
    state = score_state(x, K)
    pca_distance_error = verify_pca_parity(x, state)
    observed = h1_max(state["scores"])
    rng = np.random.default_rng(seed)
    vals = []
    for i in range(draws):
        if time.monotonic() - start > max_seconds:
            raise TimeoutError(f"Exceeded computation budget at {i}/{draws} draws")
        scores = rotated_scores(state, rng, check=(i == 0 or (i + 1) % 32 == 0))
        vals.append(h1_max(scores))
        if (i + 1) % 12 == 0:
            print(f"done: {i + 1}/{draws} covariance-preserving draws", flush=True)
    a = np.asarray(vals, dtype=float)
    sd = float(a.std(ddof=1))
    return {
        "observed_max_h1": float(observed),
        "null_max_h1": vals,
        "null_mean": float(a.mean()),
        "null_sd": sd,
        "null_min": float(a.min()),
        "null_max": float(a.max()),
        "null_exceedances": int(np.count_nonzero(a >= observed)),
        "p_one_sided_plus_one": float((1 + np.count_nonzero(a >= observed)) / (draws + 1)),
        "z_descriptive": float((observed - a.mean()) / sd) if sd > 0 else None,
        "pca_full_svd_distance_max_error": pca_distance_error,
        "num_tumor_samples": int(x.shape[0]),
        "num_features": int(x.shape[1]),
        "num_pcs": K,
        "draws": int(draws),
        "seed": int(seed),
        "seconds": round(time.monotonic() - start, 4),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", required=True, choices=["pilot", "confirm"])
    ap.add_argument("--draws", type=int)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    n = args.draws if args.draws is not None else (24 if args.mode == "pilot" else 499)
    if args.mode == "pilot" and not 2 <= n <= 32:
        ap.error("pilot permits 2-32 draws")
    if args.mode == "confirm" and n != 499:
        ap.error("confirm requires exactly 499 draws")
    seed = 20261009 if args.mode == "pilot" else 20261010
    x = load_tumor()
    result = experiment(x, n, seed, 1600 if args.mode == "pilot" else 7000)
    result.update({
        "id": "CPTAC-TUMOR-COV-NULL-001",
        "mode": args.mode,
        "status": "POST_PREREG_SENSITIVITY",
        "null": "Haar sample-centered orthogonal orbit; full feature covariance fixed",
        "assumptions": "iid matrix-normal exchangeable tumor rows, conditional on global feature selection, standardized source matrix and centered scatter; non-Gaussian/batch/patient effects not controlled",
        "input_hashes": SOURCE_HASHES,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "github_sha": os.environ.get("GITHUB_SHA", ""),
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Saved {args.out}; T_obs={result['observed_max_h1']:.6f}, "
          f"p={result['p_one_sided_plus_one']:.6f}, seconds={result['seconds']}", flush=True)


if __name__ == "__main__":
    main()
