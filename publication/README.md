# PRISM publication package

This directory freezes the publication-facing outputs for the corrected historical PRISM analysis. The authoritative inputs are under `revision_results/final_20260914/` and use the historical 125-feature FaceKit schema. Newer 120-feature FaceKit outputs are out of scope.

Run `./publication/reproduce_publication.sh` from the repository root to regenerate the numerical audit, supplementary tables, exact historical formula/dependency map, three manuscript figures, feature-correlation characterization, candidate evidence table, covariate sensitivity package, provenance record, and lint report. The command reads existing measurements and result tables; it does not rerun landmark extraction.

## Contents

- `MANUSCRIPT_RESULTS_AUDIT.md` and `manuscript_results_audit.csv`: claim-level reconstruction with source and calculation.
- `supplementary/`: all requested S1-S6 publication tables.
- `figures/`: source CSVs plus vector PDF and 350-dpi PNG outputs.
- `feature_correlation/`: the 125-by-125 reference correlation matrix, descriptive clusters, heatmap, and interpretation report.
- `feature_definitions/`: exact formulas and landmark dependencies extracted from historical FaceKit commit `e474628f3cb17a83edad59ac6e1d9e5438753032`.
- `covariate_sensitivity/`: metadata feasibility, descriptive cohort distributions, and three predefined exploratory HC3 models.
- `CANDIDATE_PHENOTYPE_EVIDENCE.csv`: evidence synthesis with literature and clinical-review status kept unresolved where appropriate.
- `REVIEWER_METHODS_AUDIT.md`: skeptical methods review and legacy-contamination inventory.
- `provenance.yaml`: revisions, source hashes, script hashes, software versions, and limitations.
- `PUBLICATION_LINT_REPORT.md`: checks of the journal manuscript and publication files.

The package treats global Benjamini-Hochberg families for pooled, matched, and interaction analyses as distinct. Photograph selection is a sensitivity analysis, not independent replication. Normalized feature ratios are not interpreted as absolute anatomical measurements.

## Deterministic CI generation

`generate_provenance.py` uses real UTC time during ordinary local runs. When `SOURCE_DATE_EPOCH` is set, it uses that timestamp instead. CI requires exact generated structure and text plus tight numeric equivalence for CSV outputs; figure files are regenerated but are not compared byte-for-byte because Matplotlib rendering and floating-point serialization can differ across operating systems. `PUBLICATION_SOURCE_COMMIT` and `PUBLICATION_PACKAGE_COMMIT` are optional explicit inputs; the package commit remains `null` until a release commit or tag can be named without creating a self-reference.

## Figure-only generation

Run `.venv/bin/python publication/scripts/generate_figures.py` to regenerate Figures 1 and 3, the new three-panel `fig2_demographic_analyses`, and two clearly labeled code-drawn schematics, without regenerating statistical analyses or tables. Use `--figures 1 2 3` for numbered figures only, or select individual figures. Numbered PDFs are synchronized to `manuscript/figures/`.

The legacy `fig2_reference_variation` files and their existing LaTeX caption are preserved. Adopting the A/B/C graphic requires changing the manuscript figure path and caption together; manuscript text is not changed by the figure generator. See [figure report and suggested captions](figures/FIGURE_GENERATION_REPORT.md) for source lineage, verified counts, output inventory, visual checks, and author integration steps.

Publication CI uses the checksum-pinned historical source in `historical_facekit/` and skips the application submodule, whose recorded commit is no longer fetchable upstream. Formula extraction verifies the preserved source against the exact historical SHA-256 before evaluation. `PUBLICATION_FACEKIT_REPOSITORY` optionally selects a Git checkout of that same revision; its source must pass the identical hash check. See the snapshot README for provenance and license.
