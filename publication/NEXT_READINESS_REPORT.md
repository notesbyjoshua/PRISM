# PRISM next publication-readiness report

## Executive status

- Numerical reproducibility: PASS. The existing manuscript audit reports 90 PASS, 0 FAIL, and 0 unresolved claims.
- Historical-method reproducibility: PASS for the 125-feature schema. Exact formulas, source locations, and landmark dependencies were recovered from FaceKit commit `e474628f3cb17a83edad59ac6e1d9e5438753032`.
- Data-sharing/privacy status: Proceeding under the user's instruction to assume no privacy concerns. The repository still needs a formal data/code-availability statement and a recorded institutional decision for submission.
- Covariate sensitivity: COMPLETE for three predefined highlighted examples using sex label, broad demographic label, yaw, pitch, and roll. Exact-age adjustment was not possible.
- Clinical literature review: INCOMPLETE. Candidate classifications and manuscript references still require verified literature and clinical review.
- Institutional approvals/statements: INCOMPLETE. Ethics, funding, data/code availability, and approved sharing language remain unresolved.
- Manuscript completeness: INCOMPLETE. The Introduction, Discussion, Conclusion, bibliography, and author-dependent statements remain placeholders.
- CI status: CONFIGURED. The workflow regenerates the package, enforces the historical schema, runs technical lint, and checks exact structure plus cross-platform numeric equivalence. A remote GitHub run is still required after push.
- Submission ready: NO

## Resolved since previous review

- The exact historical source revision is available in the FaceKit submodule and was used directly. No current 120-feature implementation was substituted.
- All 125 corrected-result features have one exact formula row, one source location, and a landmark-dependency record. No newer-only feature entered the table.
- The exact Supplementary Table S6 is separate from the earlier high-level S6 pending author review.
- Formula-derived feature dependencies and pairwise landmark/normalizer overlap are documented without changing any historical correction family.
- Common metadata feasibility and descriptive disease/reference distributions are documented.
- Three highlighted associations were analyzed in a predefined exploratory OLS family with HC3 robust standard errors and three-test Benjamini–Hochberg correction.
- The linter now distinguishes technical errors, warnings, and submission blockers. It reports 0 technical errors, 0 warnings, and 11 submission blockers.
- Publication regeneration now includes formula extraction and covariate sensitivity.
- Provenance distinguishes source and package commits and supports deterministic `SOURCE_DATE_EPOCH` timestamps in CI.
- A push/pull-request GitHub Actions workflow now checks regeneration, required outputs, exact 125-feature coverage, and semantic equivalence of regenerated CSVs while allowing platform-specific figure bytes.

## Remaining scientific blockers

- Disease status remains perfectly confounded with source dataset. Recorded-covariate adjustment cannot separate disease effects from acquisition, curation, or other source-dataset effects.
- The three adjusted models are exploratory sensitivity analyses, not replacements for the primary Mann–Whitney analysis and not independent replication.
- Candidate clinical interpretations, prior-report classifications, and any novelty language require verified literature review and clinical sign-off.
- Reference image-path uniqueness still does not establish participant uniqueness.

## Remaining human/institutional blockers

- Confirm the final data/code availability language and the public-release boundary for source-level and aggregate artifacts.
- Supply approved ethics, consent, exemption, and oversight language.
- Supply funding and grant information or an approved no-funding statement.
- Approve conflicts-of-interest language and author contributions.
- Confirm final author list, affiliations, acknowledgments, and journal-specific disclosures.

## Remaining manuscript blockers

- Complete the Introduction with verified citations.
- Complete the Discussion with clinical review and appropriately limited interpretation.
- Complete the Conclusion.
- Populate and cite `manuscript/references.bib`.
- Resolve the abstract's author-input placeholder.
- Replace all institutional and author-confirmation placeholders.
- Compile the manuscript in a LaTeX environment and review references, citations, figures, and layout. Local LaTeX compilation was not available during this pass.
- Review and approve the exact S6 before replacing the existing high-level S6.

## Recommended next actions

1. Obtain the institutional statements and author disclosures now listed as submission blockers.
2. Perform the verified literature and clinical review for the highlighted and candidate associations.
3. Review the exact feature-definition supplement with the FaceKit author, then promote it to the manuscript supplement.
4. Draft the remaining manuscript sections only after the literature and institutional inputs are available.
5. Run and review the GitHub Actions publication workflow, then compile the manuscript in the intended submission environment.
