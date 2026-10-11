# v0.3.0 — Completed matched CPTAC tumor-only null comparison

**Release date:** 2026-10-10. **Scope:** release of the fully executed and independently replayed `CPTAC-TUMOR-MATCHED-NULL-002` sensitivity study, completed after `v0.2.0`. The original manuscript, July preregistration, earlier negative results, and v0.2.0 archived files remain unchanged.

## Question and final research decision

For the **same 110 CPTAC clear-cell renal cancer tumor observations**, the same 2,000 selected and standardized protein features, identical PCA with 50 components, Euclidean metric and longest finite Vietoris–Rips one-dimensional persistence interval, is the observed topology unusual under two different null models?

| Measure | Independent-feature permutations | Covariance-preserving Haar rotations |
| --- | ---: | ---: |
| Observed maximum H1 | 5.3542976379 | 5.3542976379 |
| Null count | 499 | 499 |
| Null mean maximum H1 | 2.8410664216 | 4.4693335569 |
| Nulls as extreme as observed | 0 | 45 |
| One-sided Monte Carlo p, plus-one corrected | **0.002** | **0.092** |
| Conditional decision at p ≤ 0.05 | Reject null | Do not reject null |

**Conclusion:** The observed H1 persistence is unusually strong compared with independently permuted protein values, but **does not significantly exceed the specified covariance-preserving reference**. The signal's interpretation therefore depends on what structure the null model preserves. This is consistent with second-order covariance structure being sufficient under the Haar null; it is **not proof that covariance is the only cause** or of a distinct biological topological mechanism. The null models differ in more than one assumption. All 499 newly generated null values were independently reproduced in the same order, and a permanent source/summary/hash regression check passed.

## Exactly what is included

- [Frozen experiment design](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/audit/cptac_matched_null/PROTOCOL.md), specified before the 499 newly permuted feature-column draws, explicitly acknowledging prior awareness of the covariance-preserving result.
- [Scientific report](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/results/cptac_tumor_matched_null_20261010/REPORT.md) with results, uncertainty, assumptions, and the fact that this is an outcome-informed **post-preregistration** follow-up.
- [Full 499-draw raw record](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/results/cptac_tumor_matched_null_20261010/499_draws.json); verified SHA-256 `7f64f411c840561fb9468a7edf887836b4a9b93b5d838ad386f43211f5e87b09`.
- [Implementation](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/audit/cptac_matched_null/experiment.py), [synthetic tests](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/audit/cptac_matched_null/test_experiment.py), and [read-only evidence checker](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/audit/cptac_matched_null/verify_evidence.py).
- Original [499 covariance-preserving null draws](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/results/cptac_tumor_covariance_null_20261009/499_draws.json), unchanged.
- Updated [README](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/README.md), [limitations](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/docs/LIMITATIONS.md), and [usage instructions](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/docs/USAGE.md).
- Preserved [statistical correction notice](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md), [preregistration](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/prereg/PREREGISTRATION.md), historical manuscripts, tracked-file inventory, and SHA-256 registry.
- [Separate future-validation phase](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md) is a **plan only**, not a claim of completed validation.

## Execution and independent replay

Primary complete 499-draw execution: [GitHub Actions run 38099020457](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38099020457). Independent complete seeded replay: [run 38099306068](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38099306068). The frozen numerical array was reproduced exactly. Final merged experiment integration checks: [run 38099763640](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38099763640).

Reproduce or verify without rerunning the research:

```bash
python audit/cptac_matched_null/verify_evidence.py
```

Run nine control tests:

```bash
python -m unittest discover -s audit/cptac_matched_null -p 'test_*.py' -v
```

Those commands require the pinned Python dependencies in `requirements.txt`. The full 499-null simulation is archived; **no fresh computation is required** merely to read or cite these results.

## Known scientific boundaries

The compared null models do **not** preserve the same properties: column shuffles preserve empirical protein marginals but break dependence; Haar rotations preserve covariance but change empirical marginals and potentially subject-distribution structure. Patient/aliquot independence, batch effects, exchangeability, mixed-population panel selection and transfer to another cohort remain untested. The old class-residualization permutation p-values withdrawn in the correction notice remain invalid. No causal biological-topology discovery is asserted.

## DOI and historical preservation

- Prior immutable GitHub release `v0.2.0`: [GitHub release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0), [Zenodo DOI 10.5281/zenodo.23272317](https://doi.org/10.5281/zenodo.23272317).
- July original `v0.1.0-biorxiv`: [Zenodo DOI 10.5281/zenodo.21287944](https://doi.org/10.5281/zenodo.21287944).
- All-version Zenodo concept DOI: [10.5281/zenodo.21287943](https://doi.org/10.5281/zenodo.21287943).
- New `v0.3.0` DOI: **unknown until Zenodo processes the GitHub Release**. Use the fixed GitHub `v0.3.0` tag and commit in the interim; do not assign the old version DOI to this release.
- This is a software and numerical-evidence archive. The [LICENSE](https://github.com/smaniches/pca-topology-crossdomain-replication/blob/v0.3.0/LICENSE) grants MIT rights to the code within its stated scope, not a blanket license to manuscript text and figures.

**Closure:** Matched-null computation, replay, code, and report are finished. Any independent-cohort validation is new research requiring its own controls, study identity, and evidence record. No other experiments are silently implied by publication of this release.
