# Separate main-analysis panels

Run from the repository root:

```sh
.venv/bin/python research/main_analysis_panels/generate_facekit_panels.py --panel a
.venv/bin/python research/main_analysis_panels/generate_facekit_panels.py --panel b
```

Use `--panel both` to generate both separate figures. There is no synthetic/test-mode fallback. The starter script mentioned by the task was unavailable, so this implementation follows the supplied specification directly.

The script verifies the locked experimental primary manifest SHA256 and study ID, 9,600 unique tests, 80 eligible cells, 67 syndromes, 120 features per cell, 736 eligible patients and 982 original discoveries. It checks formulas/regions against the current feature definitions and independently reconstructs the original global BH adjustment as an audit; it never substitutes new q values into a plot. Selection uses only patient counts and category coverage, with the prescribed ties. It occurs before discovery counts are computed. The displayed 15-feature list is a constant declared before data loading.

Panel A retains all cells, using recorded-category colors and outlined syndrome-specific shapes for the five selected syndromes. Exact coincident points can overlap; the private table preserves all cells. Panel B uses the same selected five and every eligible category cell, with original signed g and global-q dots, a symmetric ±2 scale, and gray for unsupported/undefined effects. Numerical effects outside the scale are retained in the exported table.

All real-data figures, displayed values, captions and run manifests stay under verified ignored `local_review/main_analysis_panels/`. Outputs are 300-DPI PNG plus vector PDF/SVG in Times New Roman where available. No per-patient data, photographs, alternative normalization, regression, interaction tests or manuscript edits are involved. No release/commit is performed; PI/privacy review precedes any publication of small-cohort aggregates.
