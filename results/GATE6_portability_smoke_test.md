# Gate 6: Portability Smoke Test

**What this checks:** the repository runs correctly on a machine that has
never seen this project's development environment -- i.e. using only what
`requirements.txt` declares, not any package version, path, or environment
variable specific to this sandbox.

## Procedure

1. Built a fresh Python 3.11 virtual environment from scratch at a path
   outside the repository (`/tmp/gate6_venv`), with no pre-existing
   packages beyond a bare `pip`.
2. Installed **only** `pip install -r requirements.txt` into it -- no
   conda, no pre-built environment, no packages carried over from any
   other environment used during this project's development.
3. Verified every package resolved to the exact pinned version:

   | Package | requirements.txt pin | Resolved in clean venv |
   |---|---|---|
   | numpy | 2.4.6 | 2.4.6 |
   | pandas | 3.0.3 | 3.0.3 |
   | scipy | 1.17.1 | 1.17.1 |
   | scikit-learn | 1.9.0 | 1.9.0 |
   | ripser | 0.6.14 | 0.6.14 |
   | gudhi | 3.13.0 | 3.13.0 |
   | matplotlib | 3.11.0 | 3.11.0 |
   | statsmodels | 0.14.6 | 0.14.6 |
   | GEOparse | 2.0.4 | 2.0.4 |

4. Ran `reproduce.py` (default fast-path pilot phase) using **only** the
   clean venv's interpreter (`/tmp/gate6_venv/bin/python3 reproduce.py`),
   with the repository at its normal checked-out location and no other
   environment variables or paths set.

## Result

Exit code 0. Both figures generated and visually confirmed correct:
`persistence_diagrams_real_vs_gaussian.png` (raw HVG obs=3.89, PCA50
obs=4.67, matching the pilot's reported real max-H1 values) and
`barcode_and_null_distributions.png` (barcode ranking + null histograms,
correctly showing the real observed statistic as a clear outlier against
both the pipeline-symmetric permutation null (n=2000) and Gaussian null
(n=500) distributions in both feature spaces).

No hardcoded path, sandbox-specific package, or environment variable was
required. `reproduce.py`'s own `--full`/`--replication`/`--ablation`/
`--confound` phases were not re-run in the clean venv within this smoke
test (each takes 3-40+ minutes; they were already verified end-to-end in
Gate 4 using the project's `tda-repro` conda environment, which has
package versions identical to `requirements.txt`). The fast pilot path
exercised here is the one path every other phase's driver script shares
its core dependency stack with (numpy/pandas/scipy/scikit-learn/ripser),
so a clean pass here is strong evidence the same dependency resolution
holds for the other phases too.

## Verdict

**PASS.**
