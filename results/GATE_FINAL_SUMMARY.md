# Final Repo-Completeness-Gate Audit Summary

**Repository:** https://github.com/smaniches/pca-topology-crossdomain-replication
**Final commit:** `8714e62`
**Audit date:** 2026-07-09

## Gate-by-gate result

| Gate | Result | Evidence |
|---|---|---|
| 1. File Inventory | **PASS** | `git ls-files` (125 tracked files) matches on-disk contents exactly, modulo gitignored build/cache byproducts (LaTeX intermediates, `__pycache__`, fetched raw-data caches) — verified via `git check-ignore -v` on every stray file found. |
| 2. Standard Files Present | **PASS** | README.md, LICENSE, CITATION.cff, requirements.txt, reproduce.py, MANIFEST.md, checksums.sha256, prereg/ all present and current. |
| 3. No Hardcoded Absolute Paths | **PASS** | All analysis scripts use relative/argparse-configurable paths; the one prior `host.artifact_path()` sandbox dependency (GSE146889 confound script) was removed in an earlier commit by bundling its two small reused-null pickles directly into the repo. |
| 4. Reproduce Script Runs End-to-End | **PASS** | `reproduce.py` extended this session to cover all 5 phases (pilot fast/full, 3-cohort replication, 18-config ablation, 4-cohort confound audit). The pilot's full pre-registered-rigor run (n_gaussian=500, n_permutation=2000) reproduces `max_H1_persistence=3.886673` exactly and all four null-model z-scores. The other 3 phases were verified end-to-end at reduced draw count (real-data statistics deterministic and confirmed exact regardless). One transient GDC-API network failure during a concurrent test run was confirmed non-reproducible on retry (isolated re-run: exit 0). |
| 5. Number Provenance | **PASS** | `results/GATE5_number_provenance.md` traces 20+ headline manuscript numbers to source file/column. Two real issues were caught and fixed during this gate's own verification (see "Issues found and fixed" below), not merely disclosed. |
| 6. Portability Smoke Test | **PASS** | Built a venv from `requirements.txt` alone (`/tmp/gate6_venv`), ran `reproduce.py` and the fixed GSE81089 confound-audit script from it — exit 0, all reported statistics match. |
| 7. Checksum Registry | **PASS** | `checksums.sha256` (125 entries) verified self-consistent (`sha256sum -c` exit 0) and disk-vs-git diff empty except gitignored byproducts. |

## Issues found and fixed during this audit (not just disclosed)

1. **Data/Code Availability section was stale.** The manuscript's own Data and
   Code Availability paragraph still disclosed a code-packaging gap that had
   already been closed in an earlier session. Rewrote it to describe the
   actual repository state, and to precisely distinguish which phases were
   re-verified at full pre-registered draw counts (pilot only) versus reduced
   draw counts (replication/ablation/confound), rather than rounding up to
   "fully reproduced" across the board.

2. **Bibliography rendering bug.** An earlier citation-verification fix had
   placed reasoning text in a bibtex `note` field, which renders directly in
   the printed bibliography. Caught via visual proofing of the recompiled PDF;
   moved the reasoning to a `%`-prefixed bib comment (never rendered) and
   restored `note` to a normal citation tag.

3. **GATE5 provenance-table row was itself unverified.** A row in
   `GATE5_number_provenance.md` had copied an adjacent row's value with a
   `-> recheck` placeholder rather than the actual sourced value, while the
   document's own verdict claimed universal verification — an internal
   contradiction caught by independent audit. Fixed to the true sourced
   value (verified against `ablation_sweep_full_table.csv`).

4. **Real code-packaging bug: residualization-control AUC did not
   reproduce.** The manuscript's GSE81089 residualization-control claim
   (classifier AUC collapses to 0.094, permutation p=0.005) failed to
   reproduce from the packaged `confound_audit_common.py` -- it returned
   AUC=0.1108 instead. Root-caused by diffing against the original analysis's
   lineage-recovered code: the packaged `cv_auc()` function had drifted on
   two parameters simultaneously during packaging (classifier regularization
   C=1.0 vs. the original C=0.01, and AUC aggregation method). Fixed to match
   the original methodology exactly -- now reproduces AUC=0.0942. The
   associated p=0.005 permutation test had no corresponding code anywhere in
   the repository; recovered the original ad hoc test code from this
   session's own sub-agent transcript, packaged it as a new
   `auc_collapse_permutation_test()` function, wired it into the driver
   script, and verified it reproduces p_permutation=0.0050 at the correct
   200-draw count.

## What was not re-derived from scratch

The replication cohorts' (GSE146889, CPTAC-CCRCC, TCGA-LUAD), ablation
sweep's, and confound audit's specific z-scores and p-values reported in the
manuscript's Results section were established during the original full-rigor
analysis runs. This gate audit's own end-to-end re-execution ran those
phases at a reduced null-model draw count (disclosed explicitly, both in the
manuscript's Data/Code Availability section and here) -- confirming the
packaged code correctly reproduces every deterministic real-data statistic
and the general shape of the null comparisons, but not re-deriving the exact
published z/p-values at full pre-registered draw counts a second time. No
gate was silently skipped or left undisclosed.

## Final state

- All 7 gates: **PASS**
- Working tree: clean, in sync with `origin/main` at commit `8714e62`
- No open TODOs, no undefined LaTeX references, no unverifiable headline
  claims remaining in the repository as of this commit.
