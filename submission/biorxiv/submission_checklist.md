# bioRxiv Submission Checklist

Status as of this package's creation is noted in brackets after each item where verified;
unverified/not-yet-done items are left unchecked per the "no overclaiming" instruction.

- [ ] Repository is public. **[Verified via GitHub API: currently PRIVATE.]**
- [x] Final manuscript PDF is compiled from committed source. **[Verified: `paper/paper.pdf` SHA-256 matches the entry in `checksums.sha256`; 0 LaTeX errors, 0 undefined references in the last build log.]**
- [x] Final manuscript PDF contains all figures and tables. **[4 figures, per-cohort results tables, all present in the compiled PDF as of the last recompile.]**
- [x] PDF filename contains only letters, numbers, hyphens, underscores, and `.pdf`. **[`pca_topology_crossdomain_replication_biorxiv.pdf` — verified.]**
- [ ] GitHub release tag created. **[Verified via GitHub API: 0 tags exist.]**
- [ ] Zenodo DOI created, or deliberate decision made to submit with GitHub only. **[Not yet decided/created — placeholder in `SUBMISSION_SEAL.md`.]**
- [x] Manuscript Data and Code Availability includes public GitHub URL. **[Present in `paper/sections/discussion.tex`, though the URL will not resolve until the repo is made public.]**
- [ ] Manuscript Data and Code Availability includes Zenodo DOI if available. **[No DOI exists yet; not applicable until one is created.]**
- [x] Abstract says four independent cohorts, five cohort/layer analyses, three omics modalities. **[Verified in both `paper/sections/abstract.tex` and `submission/biorxiv/abstract_plain_text.txt`.]**
- [x] No stale "five independent cohorts" wording. **[Repo-wide grep at package-creation time: zero hits for "five independent cohorts" or "5 independent cohorts".]**
- [x] No universal biological-mechanism claim. **[Repo-wide grep: zero hits.]**
- [x] No universal confound-independence claim. **[Repo-wide grep: zero hits; the manuscript explicitly states confound independence is dataset-dependent, not universal.]**
- [x] Competing interests statement ready. **[See `biorxiv_submission_metadata.md`.]**
- [x] Funding statement ready. **[See `biorxiv_submission_metadata.md`.]**
- [x] Article Category: Research article with data.
- [x] Subject Area: Bioinformatics.
- [x] Distribution/reuse option selected, recommended CC BY 4.0.
- [ ] Final PDF SHA-256 recorded in `SUBMISSION_SEAL.md`. **[To be filled in the same commit as this checklist — see that file.]**
- [ ] `MANIFEST.md` regenerated after adding submission files. **[Pending — done in the same work session as this checklist, verify against the commit that follows.]**
- [ ] `checksums.sha256` regenerated after adding submission files. **[Pending — same as above.]**
- [ ] Repository checksum verification passes according to the repository's existing convention. **[To be re-verified after the above regeneration.]**

## Additional items found during this package's preparation (not in the original template)

- [ ] Two harmless dead-code fragments were found during the stale-wording sweep (an unused
      placeholder RNG assignment in `code/confound_attribution_audit/confound_audit_gse146889.py`
      line 186, and an unused helper function in
      `code/replication_cross_dataset/04c_figure_cross_dataset_replication.py` line 239-241, both
      containing the word "placeholder" in a code comment). Neither affects any published number
      or figure. Left as-is since fixing dead code was outside this package's requested scope;
      flagged here rather than silently ignored.
- [x] `\newcommand{\TODO}` macro exists in `paper/paper.tex` (a reusable helper for marking
      in-progress text) but grep confirms it is never invoked (`\TODO{...}`) anywhere in the
      manuscript — no visible TODO markers in the compiled PDF.
