# Final PRISM revision results

The approved primary analysis is complete: **1,850 unique disease patients ages 3–9**, using one eligible image per patient, compared with **942 distinct FairFace images**, all programmatically verified as age category `3-9`. Of 379 retained syndrome labels, 49 meet the unchanged main-analysis eligibility rules.

Fresh preparation reproduced the earlier corrected cohorts exactly. Matching input/code provenance justified reuse of main and Analyses A/B. Analysis C was freshly rerun for primary and alternative photo selection; it agrees with saved results. All 13 tests passed, exact historical 125-feature names and input schemas match, and independent statistical/output checks passed. No settings changed.

The primary main analysis has **6,125 tests and 811 global-BH-significant findings**, including 379 with |g|≥0.8. Available legacy tables contain 27,375 tests and 8,038 significant findings. Of those legacy significant findings, **5,577 are now ineligible**, **1,704 remain tested but are nonsignificant**, and **757 remain significant**; 54 additional common hypotheses are now significant. Effect and confidence-interval changes are supplied separately and are not inferred merely from significance loss.

Alternative photo selection changes 94 selected images and yields 834 significant main findings, with 770 shared with the primary. All-ages results remain supplementary/exploratory.

The corrected author-confirmed extraction commit, `e474628f3cb17a83edad59ac6e1d9e5438753032` (August 30, 2026), was fetched and inspected. Its historical method uses 2D roll alignment and feature-specific normalization, without yaw/pitch frontalization. See `HISTORICAL_METHODS.md` for source-grounded wording.

`FINAL_REANALYSIS_REPORT.md` explains exclusions, subgroup eligibility, A/B/C findings and provenance. Aggregate CSVs accompany complete primary and sensitivity result tables. Raw data, prior outputs and repository work were preserved. No commit, push, live-manuscript edit or comment reply was made. Patient-level manifests remain local.
