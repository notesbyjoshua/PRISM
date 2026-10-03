# PRISM publication package

This directory freezes the publication-facing outputs for the corrected historical PRISM analysis. The authoritative inputs are under `revision_results/final_20260914/` and use the historical 125-feature FaceKit schema. Newer 120-feature FaceKit outputs are out of scope.

Run `./publication/reproduce_publication.sh` from the repository root to regenerate the numerical audit, six supplementary tables, three manuscript figures, the feature-correlation characterization, candidate evidence table, provenance record, and lint report. The command reads existing measurements and result tables; it does not rerun landmark extraction.

## Contents

- `MANUSCRIPT_RESULTS_AUDIT.md` and `manuscript_results_audit.csv`: claim-level reconstruction with source and calculation.
- `supplementary/`: all requested S1-S6 publication tables.
- `figures/`: source CSVs plus vector PDF and 350-dpi PNG outputs.
- `feature_correlation/`: the 125-by-125 reference correlation matrix, descriptive clusters, heatmap, and interpretation report.
- `CANDIDATE_PHENOTYPE_EVIDENCE.csv`: evidence synthesis with literature and clinical-review status kept unresolved where appropriate.
- `REVIEWER_METHODS_AUDIT.md`: skeptical methods review and legacy-contamination inventory.
- `provenance.yaml`: revisions, source hashes, script hashes, software versions, and limitations.
- `PUBLICATION_LINT_REPORT.md`: checks of the journal manuscript and publication files.

The package treats global Benjamini-Hochberg families for pooled, matched, and interaction analyses as distinct. Photograph selection is a sensitivity analysis, not independent replication. Normalized feature ratios are not interpreted as absolute anatomical measurements.
