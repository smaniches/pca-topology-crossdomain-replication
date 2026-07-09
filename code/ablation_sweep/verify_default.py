"""
Verify the pre-registered default configuration (HVG2000_PC50_zero_np500)
against the committed ablation_sweep_full_table.csv results, to full
precision on the real-data statistics (which do not depend on permutation
count and so are unaffected by --quick runs of 02_run_sweep.py).

Expected (per repo's committed results/ablation_sweep/ablation_sweep_full_table.csv):
    real_pca_delta   = +0.7840385437011719
    z_pca_vs_pipeline = 5.037304507566495

Run: python verify_default.py
"""
import json
import os

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(_ROOT, "results")

EXPECTED_REAL_PCA_DELTA = 0.7840385437011719
EXPECTED_Z_PCA_VS_PIPELINE = 5.037304507566495
TAG = "HVG2000_PC50_zero_np500"


def main():
    results = json.load(open(os.path.join(RESULTS_DIR, "sweep_results.json")))
    row = [r for r in results if r['tag'] == TAG]
    assert row, f"config {TAG} not found in sweep_results.json -- run 02_run_sweep.py first"
    row = row[0]

    real_pca_delta = row['real_pca_delta']
    z_pca_vs_pipeline = row['z_pca_vs_pipeline']

    print(f"real_pca_delta:     computed={real_pca_delta!r}  expected={EXPECTED_REAL_PCA_DELTA!r}")
    print(f"z_pca_vs_pipeline:  computed={z_pca_vs_pipeline!r}  expected={EXPECTED_Z_PCA_VS_PIPELINE!r} "
          f"(depends on n_perm={row['n_perm_pipeline']}; only exact at n_perm=500 or full precision reruns)")

    # real_pca_delta is deterministic given the data/preprocessing/PCA seed and
    # does NOT depend on permutation/draw counts -- must match to full float precision.
    assert np.isclose(real_pca_delta, EXPECTED_REAL_PCA_DELTA, rtol=0, atol=1e-9), (
        f"real_pca_delta mismatch: {real_pca_delta} vs {EXPECTED_REAL_PCA_DELTA}")
    print("PASS: real_pca_delta matches expected value to full precision.")

    # z_pca_vs_pipeline depends on the (stochastic) permutation null; only check
    # tightly if n_perm matches the committed run's count (500), else just sanity-check sign/magnitude.
    if row['n_perm_pipeline'] == 500:
        assert np.isclose(z_pca_vs_pipeline, EXPECTED_Z_PCA_VS_PIPELINE, rtol=0.05), (
            f"z_pca_vs_pipeline drifted too far from expected: {z_pca_vs_pipeline} vs {EXPECTED_Z_PCA_VS_PIPELINE}")
        print("PASS: z_pca_vs_pipeline within 5% of expected value (n_perm=500 matches committed run).")
    else:
        assert z_pca_vs_pipeline > 3.0, "z_pca_vs_pipeline should clear the >3 sigma H0 threshold regardless of n_perm"
        print(f"NOTE: n_perm_pipeline={row['n_perm_pipeline']} != 500 (this was likely a --quick run); "
              f"z-score will differ from the full run but should still be > 3 sigma. PASS (sanity check only).")


if __name__ == "__main__":
    main()
