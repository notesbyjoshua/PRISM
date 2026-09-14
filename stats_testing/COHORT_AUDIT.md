# Cohort audit

Run from the PRISM root with `.venv/bin/python stats_testing/cohort_audit.py`.
The default output is a new timestamped directory under ignored `local_review/`.
Use `--output local_review/audit_name` to select a new directory or `--dry-run`
to calculate the audit without writing outputs. Existing output directories are
never overwritten. Run tests with `.venv/bin/python -m unittest discover -s tests -v`.

This is a descriptive audit, not a corrected statistical analysis. It does not
select a patient image, impute age, change thresholds, modify raw inputs, or
replace saved statistical results. A missing patient identifier stays unknown.
Reference file paths identify images, not people. Patient manifests and manuscript
snapshots belong in `local_review/` and must not be committed.

The audit reads both master datasets, the saved statistical configuration and
results, and the existing demographic mapping function definitions. Those functions
are isolated with the Python AST because importing the original analysis scripts
would run analyses and write outputs. A/B and C now both exclude unknown broad labels even when a subcategory is recognized.
The audit reports mapping differences to detect future divergence.

Outputs include cohort flow; source-label and age/sex distributions; per-patient
repeat-image age spans; patient metadata inconsistencies; layer-specific image and
known-patient counts; feature missingness; image-based syndrome and feature
eligibility with exclusion reasons; saved-results reconciliation; and input/code
SHA-256 hashes. `provenance.json` lists every hashed file, the current Git commits,
versions and audit assumptions. A current FaceKit commit does not establish the
version used to extract historical data.

Age summaries are conditional diagnostics: year + month/12 assumes component
fields, retains zero-zero as unresolved, and flags months outside 0–11 rather than
silently carrying them into years. Nonzero ages still need confirmation as age at
photograph. Age categories for FairFace remain categories. No age adjustment or
clinical age inference is performed.

The saved minimum-n and missingness parameters are used only to reproduce legacy
image eligibility. Nonfinite feature values are reported and excluded by the
audit; the original scripts use `dropna`, so any infinities could cause a discrepancy.
`eligibility_reconciliation.csv` identifies differences in tested syndrome-feature
pairs. The audit recognizes boolean strings explicitly and reports malformed
quality flags; this is stricter than accepting truthy strings or numeric flags.

The available master files are downstream of image export and extraction.
Upstream source-pool and rejection counts cannot be reconstructed from these
masters alone. Recover selection manifests and extraction logs before claiming a
complete source-to-analysis flow. Confirm cohort design and extraction provenance
before replacing manuscript counts, p-values, or figures.
