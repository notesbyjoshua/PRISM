# Replaced sensitivity workflow

The selected sensitivity workflow now targets the new three-category race-stratified redesign.
See [SENSITIVITY_ANALYSIS.json](../../SENSITIVITY_ANALYSIS.json) and
[the runnable methods](../../research/race_stratified_sensitivity/README.md).
`stats_testing/reference_sensitivity.py` routes to FairFace-only placebo controls, within-GMDB
syndrome-label permutations, and locked-allocation FairFace reference resampling. Outputs remain
in verified ignored `local_review/race_stratified_sensitivity/`.

## Historical pooled-study specifications

The frozen pooled primary and its published aggregate artifacts remain reproducible and unchanged.
These counts describe the older study, not the replacement sensitivity baseline:

The primary balanced reference contains 3,409 paths and yields 955 global-BH associations. The all-eligible reference (5,864 paths) yields 971, category-by-recorded-sex reference (2,114 paths) yields 945, and alternate-photo selection yields 996. Each specification retains the complete 5,760-hypothesis primary family. References and cases overlap; these are sensitivities rather than independent replications.
