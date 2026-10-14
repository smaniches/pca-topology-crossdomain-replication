"""BRCA-TRANSFER-001: prospective-within-follow-up tumor proteome null comparison.

Not a replication of a new biological mechanism or a third-party assessment.
The paired protocol was committed before this code or external-data outcome run.
Use --mode qc for provenance/eligibility only, --mode confirm for both 499-draw
reference distributions. Refuse unexpected data rather than adapting post hoc.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import helmert
from scipy.spatial.distance import pdist
from scipy.stats import beta
from sklearn.decomposition import PCA

ID = "BRCA-TRANSFER-001"
SOURCE_REV = "5863a2b27256555793b2e726d4a1cd33fbe01d82"
MIRROR = (
    "https://raw.githubusercontent.com/DominicJSykes/"
    "BreastCancerProteomics/" + SOURCE_REV + "/"
)
PROTEOME = "77_cancer_proteomes_CPTAC_itraq.csv"
CLINICAL = "clinical_data_breast_cancer.csv"
EXPECTED_BLOBS = {
    PROTEOME: "7c66480a7691c3e112d0d4bbd7c33cad11524f15",
    CLINICAL: "33fc811ca6f664f101def2ac8861bb6df1962bc4",
}
SEED_A = 2026101401
SEED_B = 2026101402
DRAWS = 499
N_FEATURES = 2000
N_PCA = 50
MIN_SAMPLES = 60
MAX_SAMPLES = 110
MAX_SOURCE_BYTES = 26_000_000
TUMOR_HEADER = re.compile(r"^([A-Z0-9]{2}-[A-Z0-9]{4})\.01TCGA$", re.I)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha(path: Path) -> str:
    return subprocess.check_output(
        ["git", "hash-object", str(path)], text=True
    ).strip()


def acquire(directory: Path) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, expected in EXPECTED_BLOBS.items():
        path = directory / name
        if not path.exists():
            request = urllib.request.Request(
                MIRROR + name,
                headers={"User-Agent": "BRCA-TRANSFER-001-evidence-fetch/1.0"},
            )
            with urllib.request.urlopen(request, timeout=90) as response:
                with path.open("wb") as f:
                    total = 0
                    while True:
                        chunk = response.read(1 << 20)
                        if not chunk:
                            break
                        total += len(chunk)
                        if total > MAX_SOURCE_BYTES:
                            raise ValueError("Downloaded source exceeds size limit")
                        f.write(chunk)
        observed = git_blob_sha(path)
        if observed != expected:
            raise ValueError(
                f"Source identity mismatch for {name}: "
                f"git-blob {observed} != {expected}"
            )
        manifest[name] = {
            "source_url": MIRROR + name,
            "git_blob_sha1": observed,
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }
    return manifest


def canonical_patient(text: object) -> str | None:
    raw = str(text).strip().upper()
    match = re.fullmatch(r"(?:TCGA-)?([A-Z0-9]{2}-[A-Z0-9]{4})(?:-[0-9A-Z]+)?", raw)
    return match.group(1) if match else None


def sample_mapping(proteome: pd.DataFrame, clinical: pd.DataFrame) -> list:
    if "Complete TCGA ID" not in clinical.columns:
        raise ValueError("Required clinical patient identifier column absent")
    clinical_ids = [canonical_patient(i) for i in clinical["Complete TCGA ID"]]
    unique_clinical_ids = {i for i in clinical_ids if i is not None}
    if not unique_clinical_ids:
        raise ValueError("Clinical source lacks recognizable TCGA IDs")
    matched = []
    for header in proteome.columns:
        match = TUMOR_HEADER.fullmatch(str(header))
        if match and match.group(1).upper() in unique_clinical_ids:
            matched.append((match.group(1).upper(), header))
    matched.sort(key=lambda x: x[0])
    ids = [i for i, _ in matched]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate patient identifiers among selected tumor columns")
    if not (MIN_SAMPLES <= len(ids) <= MAX_SAMPLES):
        raise ValueError(
            f"Expected between {MIN_SAMPLES} and {MAX_SAMPLES} "
            f"distinct matched tumor samples, found {len(ids)}"
        )
    return matched


def preprocess(proteome_path: Path, clinical_path: Path):
    proteome = pd.read_csv(proteome_path, low_memory=False)
    clinical = pd.read_csv(clinical_path, low_memory=False)
    if "RefSeq_accession_number" not in proteome.columns:
        raise ValueError("Missing RefSeq accession feature identifiers")
    mapping = sample_mapping(proteome, clinical)
    columns = [header for _, header in mapping]
    ids = [patient for patient, _ in mapping]
    refseq = proteome["RefSeq_accession_number"].fillna("").astype(str).str.strip()
    nonempty = refseq.ne("") & refseq.ne("nan")
    unique = nonempty & ~refseq.duplicated(keep="first")
    duplicates = int((nonempty & refseq.duplicated(keep="first")).sum())
    numeric = proteome.loc[unique, columns].apply(pd.to_numeric, errors="coerce")
    raw = numeric.to_numpy(dtype=np.float64)
    identifiers = refseq.loc[unique].to_numpy(dtype=str)
    measured = np.isfinite(raw)
    fractions = measured.mean(axis=1)
    elig = fractions >= 0.8
    raw, identifiers = raw[elig], identifiers[elig]
    if len(raw) < N_FEATURES:
        raise ValueError("Less than 2000 proteins meet frozen 80% finite threshold")
    medians = np.nanmedian(np.where(np.isfinite(raw), raw, np.nan), axis=1)
    if not np.isfinite(medians).all():
        raise ValueError("Non-finite protein medians after eligibility filter")
    missing_values = int((~np.isfinite(raw)).sum())
    raw = np.where(np.isfinite(raw), raw, medians[:, None])
    # Per-sample median centering first, then fixed unsupervised feature choice.
    centered = raw - np.median(raw, axis=0, keepdims=True)
    variance = centered.var(axis=1)
    good = np.isfinite(variance) & (variance > 0)
    centered, identifiers, variance = (
        centered[good], identifiers[good], variance[good]
    )
    if len(variance) < N_FEATURES:
        raise ValueError("Insufficient nonconstant proteins after normalization")
    order = np.lexsort((np.arange(len(variance)), identifiers, -variance))
    chosen = order[:N_FEATURES]
    x = centered[chosen].T.astype(np.float64, copy=True)
    x -= x.mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    if not np.isfinite(sd).all() or np.any(sd <= 0):
        raise ValueError("Non-positive selected protein sample SD")
    x /= sd
    if x.shape != (len(ids), N_FEATURES) or not np.isfinite(x).all():
        raise ValueError("Invalid normalized breast tumor matrix")
    qc = {
        "n_distinct_matched_tumors": len(ids),
        "patient_identifiers_sha256": hashlib.sha256(
            "\n".join(ids).encode("utf-8")
        ).hexdigest(),
        "tumor_identifiers": ids,
        "input_protein_rows": len(proteome),
        "duplicate_refseq_rows_excluded": duplicates,
        "eligible_at_80pct": int(elig.sum()),
        "selected_proteins": N_FEATURES,
        "eligible_protein_missing_values_median_imputed": missing_values,
        "selected_accessions_sha256": hashlib.sha256(
            "\n".join(identifiers[chosen]).encode("utf-8")
        ).hexdigest(),
        "matrix_shape": list(x.shape),
        "matrix_sha256_float64_c_order": hashlib.sha256(
            np.ascontiguousarray(x).tobytes()
        ).hexdigest(),
        "batch_covariates_available_in_this_mirror": False,
        "patient_and_batch_audit": (
            "Unique TCGA tumor columns matched against clinical IDs. "
            "Original instrument batch/channel and aliquot QC mapping NOT verified."
        ),
    }
    return x, qc


def pca_scores(x: np.ndarray, components: int = N_PCA) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 2 or not np.isfinite(x).all() or not (1 <= components < x.shape[0]):
        raise ValueError("Bad PCA input or component count")
    xc = x - x.mean(axis=0, keepdims=True)
    gram = xc @ xc.T
    values, vectors = np.linalg.eigh(gram)
    keep = np.argsort(values)[-components:][::-1]
    if np.any(values[keep] <= 1e-11):
        raise ValueError("Rank-deficient PCA leading components")
    return vectors[:, keep] * np.sqrt(values[keep])[None, :]


def assert_pca_parity(x: np.ndarray, scores: np.ndarray) -> float:
    reference = PCA(n_components=N_PCA, svd_solver="full").fit_transform(x)
    err = float(np.max(np.abs(pdist(reference) - pdist(scores))))
    if err > 1e-7:
        raise ValueError("Gram PCA != independently fitted full SVD: " + str(err))
    return err


def feature_permutation(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n, p = x.shape
    indices = np.argsort(rng.random((n, p)), axis=0, kind="stable")
    return x[indices, np.arange(p)[None, :]]


def haar_rotate(
    x: np.ndarray, rng: np.random.Generator, subspace: np.ndarray | None = None
) -> np.ndarray:
    n, _ = x.shape
    if subspace is None:
        subspace = helmert(n, full=False).T
    mu = x.mean(axis=0, keepdims=True)
    z = subspace.T @ (x - mu)
    q, r = np.linalg.qr(rng.normal(size=(n - 1, n - 1)))
    signs = np.sign(np.diag(r))
    signs[signs == 0] = 1.0
    q *= signs[None, :]
    return subspace @ (q @ z) + mu


def h1_max(scores: np.ndarray) -> float:
    from ripser import ripser
    bars = ripser(scores, maxdim=1)["dgms"][1]
    if not len(bars):
        return 0.0
    span = bars[:, 1] - bars[:, 0]
    valid = span[np.isfinite(span)]
    return float(valid.max()) if len(valid) else 0.0


def covariance_discrepancy(x: np.ndarray, rotated: np.ndarray) -> float:
    xc = x - x.mean(axis=0, keepdims=True)
    yc = rotated - rotated.mean(axis=0, keepdims=True)
    aa = xc.T @ xc
    bb = yc.T @ yc
    return float(np.linalg.norm(aa - bb) / max(np.linalg.norm(aa), 1.0))


def clopper_pearson(k: int, n: int) -> list[float]:
    low = float(beta.ppf(0.025, k, n - k + 1)) if k else 0.0
    high = float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0
    return [low, high]


def reference_report(values: np.ndarray, observed: float, seed: int, name: str):
    k = int(np.count_nonzero(values >= observed))
    return {
        "name": name, "seed": seed, "draws": len(values), "exceedances": k,
        "p_one_sided_plus_one": (k + 1) / (len(values) + 1),
        "monte_carlo_exceedance_fraction_95pct_interval":
            clopper_pearson(k, len(values)),
        "sample_mean": float(values.mean()),
        "sample_sd": float(values.std(ddof=1)),
        "min": float(values.min()), "max": float(values.max()),
        "ordered_draws": values.tolist(),
    }


def confirm(x: np.ndarray) -> dict:
    if x.shape[1] != N_FEATURES or not (MIN_SAMPLES <= x.shape[0] <= MAX_SAMPLES):
        raise ValueError("Confirmed cohort size or feature dimension invalid")
    observed_scores = pca_scores(x)
    parity_observed = assert_pca_parity(x, observed_scores)
    observed = h1_max(observed_scores)
    a_rng, b_rng = np.random.default_rng(SEED_A), np.random.default_rng(SEED_B)
    n = len(x)
    u = helmert(n, full=False).T
    a, b = np.empty(DRAWS), np.empty(DRAWS)
    parity_feature = None
    parity_haar = None
    covariance_relative_error = None
    feature_marginal_preserved = False

    for i in range(DRAWS):
        scrambled = feature_permutation(x, a_rng)
        if i == 0:
            feature_marginal_preserved = bool(np.array_equal(
                np.sort(scrambled, axis=0), np.sort(x, axis=0)
            ))
            if not feature_marginal_preserved:
                raise ValueError("Feature null changed a protein marginal")
        scores = pca_scores(scrambled)
        if i == 0:
            parity_feature = assert_pca_parity(scrambled, scores)
        a[i] = h1_max(scores)
        rotated = haar_rotate(x, b_rng, subspace=u)
        if i == 0:
            covariance_relative_error = covariance_discrepancy(x, rotated)
            if covariance_relative_error > 1e-9:
                raise ValueError("Haar null changed sample-centered covariance")
            if not np.allclose(x.mean(axis=0), rotated.mean(axis=0), atol=1e-9):
                raise ValueError("Haar null changed feature means")
        rotated_scores = pca_scores(rotated)
        if i == 0:
            parity_haar = assert_pca_parity(rotated, rotated_scores)
        b[i] = h1_max(rotated_scores)
        if not np.isfinite(a[i]) or not np.isfinite(b[i]):
            raise ValueError("Nonfinite topological statistic")
        if (i + 1) % 50 == 0:
            print(f"completed {i+1}/{DRAWS} draws in each prespecified arm", flush=True)

    arm_a = reference_report(a, observed, SEED_A, "independent-column permutation")
    arm_b = reference_report(b, observed, SEED_B, "sample-space covariance-preserving Haar")
    reject_a = arm_a["p_one_sided_plus_one"] <= 0.05
    reject_b = arm_b["p_one_sided_plus_one"] <= 0.05
    decision = (
        "same_descriptive_threshold_pattern" if reject_a and not reject_b
        else "both_reject" if reject_a and reject_b
        else "haar_only_rejects" if reject_b
        else "neither_rejects"
    )
    return {
        "observed_max_finite_h1": observed,
        "independent_feature_reference": arm_a,
        "covariance_preserving_reference": arm_b,
        "descriptive_0_05_pattern": decision,
        "checks": {
            "pca_observed_distance_max_abs_error": parity_observed,
            "pca_first_feature_null_distance_max_abs_error": parity_feature,
            "pca_first_haar_null_distance_max_abs_error": parity_haar,
            "haar_first_draw_relative_covariance_error": covariance_relative_error,
            "feature_first_draw_exact_marginal_preservation": feature_marginal_preserved,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("qc", "confirm"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="brca_transfer_") as rawdir:
        provenance = acquire(Path(rawdir))
        x, qc = preprocess(Path(rawdir) / PROTEOME, Path(rawdir) / CLINICAL)
        payload = {
            "experiment": ID, "protocol_class": "pre-outcome-specified-transfer",
            "status": "DATA_ELIGIBILITY_VALIDATED" if args.mode == "qc"
                      else "COMPUTED_PENDING_BATCH_AUDIT",
            "original_study": "PDC000173",
            "source_provenance": provenance,
            "cohort_qc": qc,
            "design": {
                "pca_components": N_PCA, "selected_protein_features": N_FEATURES,
                "reference_draws_per_arm": DRAWS,
                "h1_statistic": "longest finite Vietoris-Rips H1 on Euclidean PCA50",
            },
            "runtime": {
                "python": platform.python_version(), "numpy": np.__version__,
                "pandas": pd.__version__, "scipy": scipy.__version__,
                "platform": platform.platform(),
                "github_sha": os.environ.get("GITHUB_SHA", ""),
                "github_run_id": os.environ.get("GITHUB_RUN_ID", ""),
            },
            "limitations": [
                "Independent human cohort, but not independently replicated by an outside lab.",
                "Mirror provenance confirmed to pinned Git blob, original iTRAQ batch and aliquot QC not confirmed.",
                "Multiple tumor subtypes and instrument/batch confounding not adjusted.",
                "These two null models retain different invariants; difference of p-values does not prove covariance causation.",
                "Feature panel selected on observed tumor data before conditional reference generation.",
            ],
        }
        if args.mode == "confirm":
            payload["computed"] = confirm(x)
        payload["runtime"]["elapsed_seconds"] = round(time.monotonic() - started, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "RESULT_RECORD: " + args.mode + " " +
        json.dumps({
            "n": payload["cohort_qc"]["n_distinct_matched_tumors"],
            "input_sha256": provenance[PROTEOME]["sha256"],
            "result": payload.get("computed", {}).get("descriptive_0_05_pattern"),
            "seconds": payload["runtime"]["elapsed_seconds"],
        }), flush=True
    )


if __name__ == "__main__":
    main()
