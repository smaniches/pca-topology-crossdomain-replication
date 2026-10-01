# Derived-output reconciliation

**Author:** Santiago Maniches (ORCID: 0009-0005-6480-1987)  
**Scope:** GSE146889 clean-lineage confound-attribution audit  
**Status:** comparison-definition correction; no experiment rerun

The clean-lineage workflow correctly recomputed the raw-data, PCA, persistent-homology, null, residualization, cocycle-support, and bootstrap quantities. During review, two mismatched definitions were found only in the derived historical-comparison layer:

1. The historical report's tumor/normal CV AUC used logistic regression with `C=1.0`, but the first aggregate report compared it with the clean rerun's newly added `C=0.01` sensitivity.
2. The historical report used a two-sided Fisher exact test for cocycle-support composition, but the first aggregate report compared it with the clean rerun's one-sided enrichment test.

Both variants had already been computed from the same clean rerun. The reconciled artifacts select the matched historical definitions for reproduction claims and retain the alternatives as explicitly labelled sensitivity analyses.

No raw input, preprocessing choice, selected HVG set, PCA coordinates, persistence diagram, observed max-H1 value, Monte Carlo draw, bootstrap draw, or cocycle support was changed. The full experiment was not rerun.

After matched-definition reconciliation:

- all deterministic max-H1 values reproduce exactly;
- the historical `C=1.0` intact and residualized AUC values reproduce exactly;
- the historical two-sided Fisher p-values reproduce exactly;
- regenerated stochastic quantities preserve the same scientific decisions;
- the historical qualitative GSE146889 verdict remains sustained.

`verify_reconciliation.py` is deliberately read-only. It independently reconstructs the matched historical AUC and Fisher definitions from the committed clean-run machine-readable object, checks deterministic max-H1 values against the historical tables, verifies the bootstrap and residualization decision boundaries, and asserts that the committed reconciled comparison table and machine-readable claims agree. It does not rewrite any result.

The exact workflow used for the successful full clean-lineage execution is archived as `EXECUTION_WORKFLOW.yml` outside `.github/workflows/`; it is provenance, not an active automation. `results/SHA256SUMS` is generated only after the final reconciled artifacts are fixed and is verified by ordinary read-only CI.

## Review-driven integrity corrections

The frozen protocol is preserved as historical evidence. Its AUC annotation labelled `C=0.01` as the historical implementation, but the historical result table actually used `LogisticRegression(C=1.0)`; `C=0.01` is therefore retained only as a regularization sensitivity in the reconciled outputs.

The original preparation code also reported a permutation p-value for the class-mean-residualized AUC while holding a label-derived residualization/PCA transform fixed. That null is not exchangeable with the observed statistic because the transform was itself estimated from the true labels. The p-value is therefore removed from promoted evidence rather than replaced by an unexecuted statistic. No headline claim in this audit depends on that p-value.

The pre-reconciliation clean-run commit `979908f3c728b60add63fa26d0492d79630e2161` is the canonical source for `real_results`. A direct field-by-field comparison showed that later reconciliation edits unintentionally changed 112 raw machine-readable fields, including quartile membership and confound-spectrum values, despite the reconciliation scope stating that no experiment rerun occurred. The current repair restores the entire `real_results` object from that clean-run commit, then removes only the invalid residualized-AUC permutation p-value described above.

The retained GitHub Actions result artifact from run `32537350277` (artifact `9465987835`, archive digest `sha256:c6514242e202e0a891d40c1b5dba760ab0ff343a2fbf9c4ed0f64691b624693d`) independently preserves the original clean execution outputs. The pre-reconciliation commit confirms that `primary_results.json` was strict JSON and that q1-q4 were disjoint 44-sample quartiles partitioning indices 0..175.

`EXECUTION_WORKFLOW.yml` intentionally remains the exact archived workflow that ran, including its original active path `.github/workflows/gse146889-clean-lineage.yml`. It is provenance, not a runnable current workflow. Current `prepare.py` enforces the audited GEO SHA-256 by default, so future executions fail closed on a substituted payload.

### Original-run versus reconciled hashes

| File | Original clean-run SHA-256 | Reconciled SHA-256 |
| --- | --- | --- |
| `REPORT.md` | `512ff3f9ec10e009b787956955b71de5616c3d9f6f3d43a23d162f843f18e28b` | `b8b8227fde704d2f7bc98b2a28f87485ce6b9cdaf7fd9a0ec6288b11778a1473` |
| `historical_comparison.csv` | `94448e8b2c06b65da8f6cb3083418bd069fef892e0188e5247e5bf057b6890d7` | `d5b0fbe614f954026b78d7513631f77230d52d0434008fd344a12760ff637d54` |
| `primary_results.json` | `64d0f3e580f8c18bc1f106bf9caca42914cbe39023b704a51ce7c514ff75b2db` | `9ff196bcdbe98ef3d8d705216c10993e8c8486516cc5be7ac413497f4fec6e3a` |
| `residualization.csv` | `f05be71f5df2fb1936f53d33b84b12edb8bae00bee8b0269afb1efe4baf9084e` | `35795eff9ace14da5448acb7e79c8538e2251fd060aa2ac7eb0977e8b663eb46` |
| `run_provenance.json` | not self-hashed in the original output map | `086225182f4806268fb221b2bbad29143d91661dbd1c6813c38a2239c7a51952` |

The provenance file is intentionally preserved as the immutable record of the original clean run; the reconciliation record and integrity registries document later derived-presentation repairs explicitly.
