#!/usr/bin/env python3
"""
Single entry point to reproduce the pilot-cohort (GSE81089) figures and,
optionally, the full pilot statistics (including ~2500 persistent-homology
computations, which takes on the order of an hour on a single machine).

SCOPE (see README.md "Known gaps" for the full disclosure): this script
reproduces the pilot cohort only. The three confirmatory replication cohorts
(GSE146889, CPTAC-CCRCC, TCGA-LUAD) each have their own standalone driver
under code/replication_<COHORT>/, the 18-configuration ablation sweep is
under code/ablation_sweep/, and the confound-attribution audit (all 4
cohorts) is under code/confound_attribution_audit/ -- none of those are
wired into this entry point yet, so run them directly per their own
README.md if you need those results reproduced.

Usage:
    python3 reproduce.py                 # fast: regenerate the 2 pilot figures
                                          # from bundled checkpoints (~1 min)
    python3 reproduce.py --full          # slow: also re-run the full pilot
                                          # statistics from raw GEO data
                                          # (2000 + 500 permutation draws,
                                          # ~20-40 min depending on hardware)
"""
import os
import sys
import subprocess
import argparse

_ROOT = os.path.dirname(os.path.abspath(__file__))
_CODE = os.path.join(_ROOT, "code")
_OUT = os.path.join(_ROOT, "reproduce_output")


def run_script(script_name, cwd):
    path = os.path.join(_CODE, script_name)
    print(f"\n=== Running {script_name} (cwd={cwd}) ===")
    result = subprocess.run([sys.executable, path], cwd=cwd)
    if result.returncode != 0:
        print(f"FAILED: {script_name} exited with code {result.returncode}")
        return False
    print(f"OK: {script_name}")
    return True


def check_expected_files(cwd, filenames):
    missing = [f for f in filenames if not os.path.isfile(os.path.join(cwd, f))]
    if missing:
        print(f"FAILED: expected output file(s) not found: {missing}")
        return False
    print(f"OK: all expected output files present: {filenames}")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full", action="store_true",
        help="Also re-run the full pilot statistics from raw GEO data "
             "(slow: ~20-40 min, re-fetches GSE81089 and runs 2000+500 "
             "permutation/Gaussian-null persistent-homology computations).",
    )
    args = parser.parse_args()

    os.makedirs(_OUT, exist_ok=True)
    ok = True

    # Fast path: regenerate the 2 pilot figures from bundled checkpoints
    # (results/pilot_GSE81089/checkpoints/*.pkl).
    ok &= run_script("04a_figure_persistence_diagrams.py", cwd=_OUT)
    ok &= check_expected_files(_OUT, ["persistence_diagrams_real_vs_gaussian.png"])

    ok &= run_script("04b_figure_barcode_null_distributions.py", cwd=_OUT)
    ok &= check_expected_files(_OUT, ["barcode_and_null_distributions.png"])

    if args.full:
        # Slow path: re-fetch GSE81089 from GEO and recompute all statistics
        # from scratch (does not depend on the bundled checkpoints).
        ok &= run_script("03_statistics_and_results_table.py", cwd=_OUT)
        ok &= check_expected_files(_OUT, ["final_results_table.csv"])
        if ok:
            import pandas as pd
            df = pd.read_csv(os.path.join(_OUT, "final_results_table.csv"))
            # Spot-check against the pre-registered, previously-reported numbers
            # (see results/pilot_GSE81089/final_results_table.csv).
            row = df[df["comparison"] == "Real HVG-raw vs pipeline-null"].iloc[0]
            expected_max_pers = 3.886673
            actual = row["observed_or_diff"]
            if abs(actual - expected_max_pers) > 1e-3:
                print(f"FAILED: reproduced max_H1_persistence ({actual}) does not "
                      f"match the expected value ({expected_max_pers}) to tolerance")
                ok = False
            else:
                print(f"OK: reproduced max_H1_persistence ({actual}) matches "
                      f"expected value ({expected_max_pers})")
    else:
        print("\n(--full not passed: skipping the ~20-40 min full statistics "
              "re-run from raw GEO data. Pass --full to include it.)")

    print(f"\nOutputs written to: {_OUT}")
    print("REPRODUCE: " + ("SUCCESS" if ok else "FAILURE"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
