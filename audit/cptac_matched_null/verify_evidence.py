"""Read-only checks for the archived matched CPTAC tumor null result.

Never computes new null draws, changes archived evidence or infers biological
causation. The reference family is conditional on its own transformation.
"""
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NEW = ROOT / "results/cptac_tumor_matched_null_20261010/499_draws.json"
HAAR = ROOT / "results/cptac_tumor_covariance_null_20261009/499_draws.json"
EXPECTED_NEW_SHA = "7f64f411c840561fb9468a7edf887836b4a9b93b5d838ad386f43211f5e87b09"
EXPECTED_HAAR_SHA = "e45385676be011c6f2df7f6dbbe30202d5a28418bc6e490e8af176e77b89fcfa"
EXPECTED_VECTOR_SHA = "b3b93df2c69eb62b95964c422dd6b88886cd99df8b5398d2b27fa4a67b945cb4"
OBS = 5.354297637939453


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_locked(path, expected_digest):
    raw = path.read_bytes()
    actual = digest(raw)
    if actual != expected_digest:
        raise ValueError(f"Evidence SHA-256 mismatch: {path}: {actual}")
    return json.loads(raw)


def check_summary(values, observed, reported, counts_key):
    if len(values) != 499 or any(not isinstance(x, (int, float)) or
                                 not math.isfinite(x) or x < 0 for x in values):
        raise ValueError("Null array is not 499 finite positive numbers")
    exceed = sum(value >= observed for value in values)
    p = (1 + exceed) / 500
    checks = [
        (reported["exceedances" if counts_key == "new" else "null_exceedances"], exceed),
        (reported["p_one_sided_plus_one"], p),
        (reported["null_mean"], statistics.mean(values)),
        (reported["null_sd"], statistics.stdev(values)),
        (reported["null_max"], max(values)),
        (reported["null_min"], min(values)),
    ]
    for actual, expected in checks:
        if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Reported summary inconsistent with unrounded values")
    return exceed, p


def main():
    new = read_locked(NEW, EXPECTED_NEW_SHA)
    prior = read_locked(HAAR, EXPECTED_HAAR_SHA)
    if new["experiment_id"] != "CPTAC-TUMOR-MATCHED-NULL-002":
        raise ValueError("Wrong experimental identity")
    if (new["sample_count"], new["features"], new["pca_components"]) != (110, 2000, 50):
        raise ValueError("Observed shape or analysis settings were altered")
    if abs(new["observed_max_h1"] - OBS) > 1e-12:
        raise ValueError("Changed observed H1")
    if abs(prior["observed_max_h1"] - OBS) > 1e-12:
        raise ValueError("Comparator uses different observed H1")
    a = new["new_reference"]
    b = new["archived_covariance_preserving_reference"]
    if a["draws"] != 499 or a["seed"] != 20261011 or b["draws"] != 499 or b["seed"] != 20261010:
        raise ValueError("One of the frozen seeds/draw counts changed")
    feature_values = a["null_max_h1"]
    if digest(json.dumps(feature_values, separators=(",", ":")).encode()) != EXPECTED_VECTOR_SHA:
        raise ValueError("Ordered source null vector SHA-256 mismatch")
    feature_exceed, feature_p = check_summary(feature_values, OBS, a, "new")
    prior_exceed, prior_p = check_summary(prior["null_max_h1"], OBS, prior, "prior")
    if feature_exceed != 0 or feature_p != 0.002:
        raise ValueError("New feature-null decision changed")
    if prior_exceed != 45 or prior_p != 0.092:
        raise ValueError("Prior Haar decision changed")
    if (b["json_sha256"] != EXPECTED_HAAR_SHA or
            b["exceedances"] != prior_exceed or
            b["p_one_sided_plus_one"] != prior_p or
            not math.isclose(b["null_mean"], prior["null_mean"], abs_tol=1e-10)):
        raise ValueError("Matched comparison does not reflect actual archived Haar source")
    if (not math.isfinite(new["checks"]["observed_metric_full_svd_max_abs_error"]) or
            new["checks"]["observed_metric_full_svd_max_abs_error"] > 1e-8 or
            new["checks"]["first_null_metric_full_svd_max_abs_error"] > 1e-8):
        raise ValueError("PCA/SVD numerical parity gate failed")
    print("MATCHED-NULL ARCHIVE: PASS")
    print("Same observed max-H1:", OBS)
    print("Feature-shuffled: exceed 0/499, p=0.002")
    print("Covariance-preserving: exceed 45/499, p=0.092")
    print("This verifies archived computations, not biological causation.")


if __name__ == "__main__":
    main()
