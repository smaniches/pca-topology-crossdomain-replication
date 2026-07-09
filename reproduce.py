#!/usr/bin/env python3
"""
Single entry point to reproduce every phase of the TOPOLOGICA PCA-topology
cross-domain replication project: the pilot cohort (GSE81089), the 3
confirmatory replication cohorts (GSE146889, CPTAC-CCRCC, TCGA-LUAD RNA-seq
+ methylation), the 18-configuration ablation sweep, and the 4-cohort
confound-attribution audit.

Each phase's underlying code lives in its own subdirectory under code/ and
can also be run standalone -- see each subdirectory's README.md for full
detail (exact expected values, tolerances, and any documented pipeline
deviations). This script is a thin, phase-selectable orchestrator over
those standalone scripts; it does not duplicate their logic.

DEFAULT (no flags): fast path only -- regenerates the 2 pilot figures from
bundled checkpoints (~1 min), the cheapest possible sanity check that the
repo is wired together correctly. Everything else is opt-in via flags,
because most phases involve re-fetching public data (tens of MB to ~150 MB)
and/or non-trivial persistent-homology compute time.

Usage:
    python3 reproduce.py                     # fast: pilot figures only (~1 min)
    python3 reproduce.py --full              # slow: full pilot statistics from
                                              # raw GEO data (2000+500 draws, ~20-40 min)
    python3 reproduce.py --replication        # re-fetch + verify all 3 confirmatory
                                              # replication cohorts (reduced draw
                                              # count by default; see --n-draws)
    python3 reproduce.py --ablation          # run the 18-config ablation sweep
                                              # (reduced draw count by default)
    python3 reproduce.py --confound          # run the 4-cohort confound-attribution
                                              # audit (reduced draw count by default)
    python3 reproduce.py --all               # everything above
    python3 reproduce.py --all --n-draws 500 # everything, at (closer to) the
                                              # pre-registered draw counts --
                                              # can take several hours total

Every phase beyond the fast default prints its own real-data statistics and
asserts them against the values in the corresponding results/*/*.md report
(see each script's docstring for exact tolerances). A non-zero exit code
from this script means at least one phase's own assertions failed --
grep the printed output for "FAILED" to find which.
"""
import os
import sys
import subprocess
import argparse

_ROOT = os.path.dirname(os.path.abspath(__file__))
_CODE = os.path.join(_ROOT, "code")
_OUT = os.path.join(_ROOT, "reproduce_output")

# Default reduced draw count for the opt-in phases (keeps --all runnable in
# minutes rather than hours; real-data statistics do not depend on this and
# are always asserted at full tolerance regardless of the count used).
_DEFAULT_QUICK_DRAWS = 20


def run(cmd, cwd, label):
    print(f"\n=== {label} (cwd={cwd}) ===")
    print("    $", " ".join(cmd))
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"FAILED: {label} exited with code {result.returncode}")
        return False
    print(f"OK: {label}")
    return True


def run_script(script_name, cwd, args=None):
    path = os.path.join(_CODE, script_name)
    return run([sys.executable, path] + (args or []), cwd=cwd, label=script_name)


def check_expected_files(cwd, filenames):
    missing = [f for f in filenames if not os.path.isfile(os.path.join(cwd, f))]
    if missing:
        print(f"FAILED: expected output file(s) not found: {missing}")
        return False
    print(f"OK: all expected output files present: {filenames}")
    return True


def phase_pilot_fast(ok):
    ok &= run_script("04a_figure_persistence_diagrams.py", cwd=_OUT)
    ok &= check_expected_files(_OUT, ["persistence_diagrams_real_vs_gaussian.png"])
    ok &= run_script("04b_figure_barcode_null_distributions.py", cwd=_OUT)
    ok &= check_expected_files(_OUT, ["barcode_and_null_distributions.png"])
    return ok


def phase_pilot_full(ok):
    ok &= run_script("03_statistics_and_results_table.py", cwd=_OUT)
    ok &= check_expected_files(_OUT, ["final_results_table.csv"])
    if ok:
        import pandas as pd
        df = pd.read_csv(os.path.join(_OUT, "final_results_table.csv"))
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
    return ok


def phase_replication(ok, n_draws):
    da = ["--n-gauss", str(n_draws), "--n-perm", str(n_draws)]
    ok &= run_script(os.path.join("replication_GSE146889", "01_fetch_preprocess_and_results_table.py"),
                      cwd=_CODE, args=da)
    ok &= run_script(os.path.join("replication_CPTAC_CCRCC", "01_fetch_preprocess_and_results_table.py"),
                      cwd=_CODE, args=da)
    ok &= run_script(os.path.join("replication_TCGA_LUAD", "01_rnaseq_fetch_preprocess_and_results_table.py"),
                      cwd=_CODE, args=da)
    ok &= run_script(os.path.join("replication_TCGA_LUAD", "02_methylation_fetch_preprocess_and_results_table.py"),
                      cwd=_CODE, args=da)
    return ok


def phase_ablation(ok, n_draws):
    ab_dir = os.path.join(_CODE, "ablation_sweep")
    ok &= run_script(os.path.join("ablation_sweep", "00_download_and_preprocess.py"), cwd=_CODE)
    ok &= run_script(os.path.join("ablation_sweep", "01_make_configs.py"), cwd=_CODE)
    ok &= run(
        [sys.executable, os.path.join(_CODE, "ablation_sweep", "02_run_sweep.py"), "--quick", str(n_draws)],
        cwd=_CODE, label="ablation_sweep/02_run_sweep.py (all 18 configs)",
    )
    ok &= run_script(os.path.join("ablation_sweep", "03_aggregate_results.py"), cwd=_CODE)
    ok &= run_script(os.path.join("ablation_sweep", "verify_default.py"), cwd=_CODE)
    return ok


def phase_confound(ok, n_draws):
    qa = ["--quick", str(n_draws)]
    ok &= run_script(os.path.join("confound_attribution_audit", "confound_audit_gse81089.py"), cwd=_CODE, args=qa)
    ok &= run_script(os.path.join("confound_attribution_audit", "confound_audit_cptac_ccrcc.py"), cwd=_CODE, args=qa)
    ok &= run_script(os.path.join("confound_attribution_audit", "confound_audit_tcga_luad.py"), cwd=_CODE, args=qa)
    ok &= run_script(os.path.join("confound_attribution_audit", "confound_audit_gse146889.py"), cwd=_CODE, args=qa)
    return ok


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true",
                         help="Also re-run the full PILOT statistics from raw GEO data "
                              "(slow: ~20-40 min, 2000+500 draws).")
    parser.add_argument("--replication", action="store_true",
                         help="Re-fetch and verify all 3 confirmatory replication cohorts "
                              "(GSE146889, CPTAC-CCRCC, TCGA-LUAD RNA-seq + methylation).")
    parser.add_argument("--ablation", action="store_true",
                         help="Run the 18-configuration ablation sweep on GSE81089.")
    parser.add_argument("--confound", action="store_true",
                         help="Run the 4-cohort confound-attribution audit.")
    parser.add_argument("--all", action="store_true",
                         help="Shorthand for --full --replication --ablation --confound.")
    parser.add_argument("--n-draws", type=int, default=_DEFAULT_QUICK_DRAWS,
                         help=f"Null-model draw count for --replication/--ablation/--confound "
                              f"(default {_DEFAULT_QUICK_DRAWS}, a fast smoke test; the "
                              f"pre-registered counts are 500-2000 and take substantially "
                              f"longer). Real-data statistics are always verified at full "
                              f"tolerance regardless of this value; only null-model z-scores/"
                              f"p-values become noisier at low draw counts.")
    args = parser.parse_args()

    if args.all:
        args.full = args.replication = args.ablation = args.confound = True

    os.makedirs(_OUT, exist_ok=True)
    ok = True

    print("### Phase: pilot cohort (GSE81089), fast path (bundled checkpoints) ###")
    ok = phase_pilot_fast(ok)

    if args.full:
        print("\n### Phase: pilot cohort (GSE81089), full statistics from raw GEO data ###")
        ok = phase_pilot_full(ok)
    else:
        print("\n(--full not passed: skipping the pilot's ~20-40 min full statistics "
              "re-run from raw GEO data.)")

    if args.replication:
        print(f"\n### Phase: 3 confirmatory replication cohorts (n_draws={args.n_draws}) ###")
        ok = phase_replication(ok, args.n_draws)
    else:
        print("\n(--replication not passed: skipping the 3 confirmatory replication cohorts. "
              "Each has its own standalone script under code/replication_<COHORT>/.)")

    if args.ablation:
        print(f"\n### Phase: 18-configuration ablation sweep (n_draws={args.n_draws}) ###")
        ok = phase_ablation(ok, args.n_draws)
    else:
        print("\n(--ablation not passed: skipping the 18-config ablation sweep. "
              "See code/ablation_sweep/README.md to run it standalone.)")

    if args.confound:
        print(f"\n### Phase: 4-cohort confound-attribution audit (n_draws={args.n_draws}) ###")
        ok = phase_confound(ok, args.n_draws)
    else:
        print("\n(--confound not passed: skipping the 4-cohort confound-attribution audit. "
              "See code/confound_attribution_audit/README.md to run it standalone.)")

    print(f"\nOutputs written to: {_OUT} (pilot fast/full path) and under each phase's own "
          f"code/<phase>/ subdirectory (replication/ablation/confound phases -- these write "
          f"their outputs alongside their own scripts, not into {_OUT}; see each phase's "
          f"printed \"Saved results to ...\" / \"Wrote ...\" lines above).")
    print("REPRODUCE: " + ("SUCCESS" if ok else "FAILURE"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
