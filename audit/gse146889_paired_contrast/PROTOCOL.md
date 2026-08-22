# GSE146889 paired-contrast topology falsification protocol

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)  
**Status:** frozen before execution  
**Purpose:** test whether the GSE146889 topological signal survives a patient-paired tumor-minus-normal representation that removes each patient's shared baseline and respects the cohort's matched sampling structure.

This is a new, prospective robustness test of the already completed clean-lineage GSE146889 audit. It is not a retroactive preregistration of the original study, is not external human validation, and does not alter the previously reported results regardless of outcome.

## 1. Scientific question

For the 85 patients with both tumor and adjacent-normal RNA-seq samples, does the point cloud of within-patient expression contrasts contain maximum finite H1 persistence larger than expected when tumor/normal orientation is randomized within each patient?

For complete pair i, define

    d_i = x_i,tumor - x_i,normal.

The primary null randomizes only the orientation of these fixed patient-level contrasts:

    d_i* = s_i d_i,  s_i in {-1,+1}, independently across pairs.

This is the matched-pair label-swap/random-sign null. It preserves each patient's pairwise contrast magnitude and all feature dependence within a pair while destroying systematic tumor-versus-normal orientation. Its randomization interpretation relies on within-pair tumor/normal label exchangeability under the sharp null; that assumption is part of the test definition, not an established biological fact.

## 2. Frozen input lineage

Primary input is the exact prepared matrix produced by the successful clean-lineage audit:

- clean-lineage workflow run: `32532718967`
- prepared artifact name: `gse146889-prepared`
- prepared artifact ZIP SHA-256: `81ffa1b920b8978cc415a189b2c2b8863c22dcf69af3a571eca630c280429c91`
- `prepared.npz` SHA-256: `75a6d4949482e4d6a06c46bbfd1b925c74e3981ae195b1b8b545af91765c2b21`
- originating raw GEO matrix SHA-256: `16e22f1285cb3960cc30151e67f8bf1ccc94a5d2e42051de9cc47ca1bbba3fa4`
- mixed feature matrix: `mixed_std`, shape 176 x 2000
- sample identifiers: `sample_names`, length 176

The 2,000-feature surface is therefore frozen before this analysis. No features are selected, dropped, reordered, or tuned using the paired outcome. It was selected unsupervised on the full mixed cohort (including the six tumor-only samples), so the paired test is conditional on this pre-existing representation rather than a pair-native feature-selection procedure.

For future reproduction after the workflow artifact expires, the same `prepared.npz` must be regenerated from the raw GEO matrix using the committed clean-lineage `prepare.py`, its frozen protocol, and the pinned environment. Execution must abort if the regenerated SHA-256 differs from the value above.

## 3. Pair construction

Sample identifiers must match exactly:

    <group>_(tumor|normal)_<integer>

The pair key is `<group>_<integer>`. A key enters the primary analysis only if exactly one tumor and exactly one normal sample are present. Ambiguous or duplicate keys abort execution.

Feasibility inspection before outcome computation established the structural expectation:

- 91 patient keys total;
- 85 complete tumor/normal pairs;
- 6 tumor-only keys;
- 0 normal-only keys.

The six incomplete keys are excluded by this rule before any paired topological statistic is computed. Their identities are structurally fixed as `MSI_MLH1G_5`, `MSI_MSH2_3`, `MSI_MSH2_8`, `MSI_MSH6_9`, `MSI_PMS2_8`, and `MSS_3`; any mismatch aborts execution. The complete-pair count must remain exactly 85.

## 4. Primary representation and statistic

1. For every complete pair, compute `d_i = mixed_std[tumor] - mixed_std[normal]` on the frozen 2,000-feature surface.
2. Standardize each of the 2,000 contrast features across the 85 pair rows with `StandardScaler` fit only to the paired-contrast matrix.
3. Fit `PCA(n_components=50, random_state=42)` under the pinned scikit-learn version.
4. Compute Vietoris-Rips persistent homology using `ripser(..., maxdim=1, coeff=2, metric="euclidean")` on the 50-dimensional PCA point cloud.
5. Primary statistic `T_obs` is the maximum finite H1 death-minus-birth persistence. If no finite H1 class exists, `T_obs = 0`.

No outcome-dependent dimensionality choice, filtration change, metric change, feature reselection, or alternative statistic is permitted for the primary test.

## 5. Primary paired sign-flip null

Generate exactly 2,000 Monte Carlo paired-randomization draws. This samples the conditional label-swap null; it does not enumerate all `2^85` possible sign configurations.

For draw j in 0..1999:

1. Generate 85 independent Rademacher signs using `numpy.random.SeedSequence([14689401, j])` and `default_rng`.
2. Multiply each complete-pair contrast row by its sign.
3. Re-fit `StandardScaler` across rows.
4. Re-fit PCA(50), random_state=42.
5. Recompute maximum finite H1 persistence.

All 2,000 draw values are preserved.

Report:

- null mean and sample SD (`ddof=1`);
- z = `(T_obs - mean_null) / sd_null`;
- one-sided Monte Carlo p = `(1 + # {T_null >= T_obs}) / 2001`;
- null 95th and 99th percentiles.

The Monte Carlo p-value is the primary inferential quantity. The z-score is a standardized effect-size diagnostic retained for continuity with the existing program; normality of the null is not assumed.

## 6. Pair-block bootstrap stability

Generate exactly 2,000 bootstrap draws.

For draw j in 0..1999:

1. Use `SeedSequence([14689402, j])` and sample 85 complete pair rows with replacement.
2. Re-fit `StandardScaler`.
3. Re-fit PCA(50), random_state=42.
4. Recompute maximum finite H1 persistence.

Report:

- percentile 95% interval for bootstrap `T`;
- percentile 95% interval for `T_boot - mean_null`, holding the sign-flip null mean fixed;
- the fraction of bootstrap draws with `T_boot > mean_null`.

The fixed-null delta interval is explicitly a stability approximation, not a joint bootstrap of observed and null processes.

## 7. Frozen decision rule

The paired audit is classified before results are examined:

### Strong paired survival

All three conditions hold:

1. one-sided sign-flip p <= 0.01;
2. z > 3;
3. the 95% pair-bootstrap interval for `T_boot - mean_null` has lower bound > 0.

### Partial / unstable paired evidence

The sign-flip p <= 0.05, but one or both of the z > 3 and bootstrap-stability criteria fail.

### No paired survival evidence

The sign-flip p > 0.05.

These labels govern only this paired robustness question. A failure cannot invalidate the already reproduced aggregate computation; it instead narrows the biological/geometric interpretation by showing dependence on between-patient structure, orientation, or paired-sample aggregation.

No threshold will be changed after seeing results.

## 8. Required outputs

Before interpretation, execution must write:

- `pair_manifest.csv`: pair key, tumor sample, normal sample, original row indices;
- `excluded_incomplete_pairs.csv`;
- `primary_result.json`;
- `sign_flip_draws.csv` with all 2,000 values;
- `pair_bootstrap_draws.csv` with all 2,000 values;
- `summary.csv`;
- `run_provenance.json` including every input/script/environment hash and seed rule;
- `paired_contrast_audit.svg` generated only from machine-readable outputs;
- `REPORT.md` with the frozen classification and explicit limitations;
- `SHA256SUMS` for every result file except itself.

The implementation must be deterministic by draw index so sharding/parallelism cannot alter the random draws. Each draw file must record its integer draw index; aggregation must assert the exact index set `0..1999` once and only once for each stochastic procedure.

## 9. Interpretation boundary

Even strong paired survival would support only the following statement:

> On the frozen GSE146889 2,000-HVG representation, patient-level tumor-minus-adjacent-normal contrasts retain H1 structure that is unusually persistent relative to within-pair label-orientation randomization and is stable under pair resampling.

It would not establish:

- a causal tumor mechanism;
- tissue-independent biology;
- a canonical biological loop;
- diagnostic or clinical validity;
- independence from every latent confound;
- replication in another cohort.

Because the 2,000-HVG surface was selected unsupervised on the original mixed cohort, a positive result would still inherit that representation choice. A future pair-native feature-selection analysis, if scientifically justified, must be specified as a separate experiment rather than selected after this result is known.

## 10. Execution boundary

This protocol must be committed on the paired-audit branch before the implementation is run or the observed paired max-H1 statistic is computed. Feasibility inspection was limited to input shapes, sample-name pairing, and package/runtime availability; no paired persistence outcome has been calculated before this freeze. The implementation may then be added in a subsequent commit. Results must be committed only to the paired-audit branch or attached as workflow artifacts. No merge, release, manuscript rewrite, or claim expansion is authorized by execution alone.
