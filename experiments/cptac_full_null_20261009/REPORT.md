# CPTAC-CCRCC full-count re-execution: verified result

**Experiment:** `CPTAC-FULLNULL-20261009`  
**Status:** completed and verified within the scope below; historical manuscript unchanged.  
**Run:** [GitHub Actions #37935680700](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/37935680700), successful on the first attempt, 2026-10-09.  
**Source:** original repository `main` at `70f7f6db91f59839eadc66bcbc1151e03e248edc`, experiment PR [#6](https://github.com/smaniches/pca-topology-crossdomain-replication/pull/6), logged workflow checkout `45f2aa0acee1262f3ae3e9cb4849d719fa16550d` (GitHub pull-request merge ref).  
**Research category:** re-execution of an already reported preregistered analysis, **not** a new independent biological cohort, prospective test, or alternative scientific implementation.

## What `n` denotes

- **Biological cohort:** n=194 actual proteomic samples (110 primary tumors, 84 normal tissues), PDC000127.
- **Feature dimensions:** 9,591 measured protein columns in the cached matrix; top 2,000 high-variance features selected and PCA(50) applied.
- **Monte Carlo draws:** N=500 Gaussian-null samples and N=2,000 pipeline-symmetric feature permutations. These are **2,500 simulated null realizations, not 2,500 new patients**.
- **Evaluation:** one real-data max-H1 statistic in raw HVG space and one in PCA50 space; each compared with both null distributions. Four statistical comparison rows, but two actual null families.

## Execution and preflight

The input matrix and labels passed the previously committed SHA-256 checks. The workflow independently checked that 194 unique matrix sample IDs exactly matched label order and that class counts were 110/84. A 2+2-draw smoke test passed real-data H1 and five-fold CV-AUC checks. The full run then completed 500/500 Gaussian and 2,000/2,000 permutation draws with seed 42, one standard GitHub-hosted Ubuntu runner and Python 3.11.17. Original code, original preprocessing, and the locked hypotheses were unchanged. The full run used `--skip-assert` after the preflight assertions; the post-run schema and numerical-invariant checks passed.

### Full-run statistics

| Comparison | Observed max-H1 | Null mean | Null SD | z | Monte Carlo p | Null draws |
|---|---:|---:|---:|---:|---:|---:|
| Raw HVG vs Gaussian | 3.832947 | 1.026919 | 0.117994 | 23.781 | 0.001996 | 500 |
| PCA50 vs Gaussian | 4.683274 | 3.111337 | 0.366643 | 4.287 | 0.001996 | 500 |
| Raw HVG vs permutation | 3.832947 | 0.985232 | 0.131341 | 21.682 | 0.000500 | 2,000 |
| PCA50 vs permutation | 4.683274 | 2.998558 | 0.356531 | 4.725 | 0.000500 | 2,000 |

**Pre-registered H0 effect-size gate:** all four z scores exceed 3.0, as previously reported.

**Pre-registered H1 direction:** real PCA-driven max-H1 increase = **0.850327**; the corresponding mean null increases are **2.084418** (Gaussian) and **2.013327** (permutation). Both exceed the real increase, sustaining the original direction.

**Historical comparison:** deterministic observed max-H1 values match the archived report to floating-point precision. Null means, SDs and z scores differ slightly. The Gaussian/PCA Monte Carlo p estimate changed from 0.003992 to 0.001996, while remaining in the same qualitative direction. No historical output was overwritten. See `summary.json` for every paired historical difference at full precision.

## Verification and artifact provenance

- Ordinary repository CI passed for the research PR.
- Experiment workflow, dataset SHA-256 checks, class-label alignment, preflight, full simulation, result-schema invariants and artifact upload all passed.
- The archived `results.csv`, `summary.json` and `pip-freeze.txt` are **byte-identical** to their original GitHub Actions artifact copies; checksums listed below.
- The workflow that performed the computation is archived as `EXECUTION_WORKFLOW.yml` outside `.github/workflows/`, so merging this report does not leave the one-time expensive experiment active.

| Archived file | Original SHA-256 |
|---|---|
| `results.csv` | `0ce1a190b90c6a33a8aa2f3c70d0b71ff88fb18feabe32da0d384ec12c8d65c3` |
| `summary.json` | `99400e216b4552216d0021e67f3fc973212df9b43cfe07b0ca4463c38c7d2115` |
| `pip-freeze.txt` | `85358288a37175713b7a87b39b02a5002183e611c4e6c3b3c3b505b99913f96b` |

The original, complete 9-file GitHub Actions artifact is `cptac-full-null-20261009-37935680700` (ID `11620714289`, ZIP digest `sha256:e813f5e05da0f0fabb7834b7bcf3fe8d4ec427733129ed49ca2fc9ee64e5463e`). [Open original run and artifact](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/37935680700). The artifact is subject to GitHub retention (original expiry 2026-11-08); the three primary machine-readable evidence records above are committed permanently once PR #6 merges.

## Limits and decisions

- This is a **same-code, same-input re-execution**, not a fully independent replication or evidence of new biological generalization.
- The independent TwoNN gate and metric-cascade selection were **not** rerun; the archived report supplies that context.
- The pre-registered 20-test Benjamini-Hochberg correction was **not** recomputed as a family-wide post-run result; do not infer a fresh corrected p-value from this report.
- The 1.000 tumor-vs-normal classifier CV-AUC is a strong label-separation signal, not evidence for biological causal interpretation of individual loops.
- Scientific state: **reproducibility strengthened; no changed scientific verdict.**
- Software state: **no change to the original computational method**; permanent provenance added, one-time workflow retired. No PyPI release applies.

The next experiment, if any, must reduce a different open scientific uncertainty rather than rerun these exact 2,500 draws again.
