# Current PRISM data layout

The validated 120-feature study reads private cohorts through `ACTIVE_ANALYSIS.json`.
Current GMDB, expanded FairFace and separate synthetic PDIDB measurements, image
inputs, checkpoints and linkage files stay in ignored private storage. They are
not copied into this public directory.

`syndrome_folder_map.csv` is retained because the current workflow uses this
condition-label taxonomy. It contains syndrome labels and condition identifiers,
not per-person observations or feature measurements.

The obsolete measurement CSVs, their row-level metadata, legacy formatting scripts
and comparison-only exports have been removed from working locations and preserved
in a verified private read-only archive. Shared GMDB Arrow inputs and the validated
120-feature measurements are retained. Existing Git objects and remote copies may
still contain retired files; see `publication/PI_REVIEW.md` for history remediation.
