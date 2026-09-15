# Corrected PRISM reanalysis

14 September 2026. Primary cohort and age handling are now author-approved. This package supersedes the provisional scope and current-version normalization wording in earlier local reports; it does not edit the manuscript.

## What changed and what was reused

The primary cohort uses supplied GMDB metadata with 3 ≤ year + month/12 < 10, the existing integer-component validation, month 0–11, and explicit exclusion of missing, invalid and zero-zero ages. No exceptional values were reinterpreted and no manual age-record review was initiated. One frontal-pass eligible image is selected per known disease patient using the existing lowest normalized absolute pose score, then lexical image ID. The alternative uses lexical image ID after the same eligibility rules. All-ages results are supplementary/exploratory only.

Fresh cohort preparation reproduced all prior cohort CSVs byte-for-byte. All 11 saved code hashes and 13 saved input hashes matched. Main statistics and A/B results were therefore reused into this new directory, as authorized, rather than described as new executions. C was freshly executed for primary and alternative cohorts to eliminate the earlier execution-time code ambiguity; results match the saved tables to numerical tolerance. No statistical settings changed. Raw data and original result directories remain intact.

The original FairFace reference predates importer fixes. Fixing its importer did not repair historical files. The existing derivative removes 58 identical duplicate-path rows from 1,000, yielding 942 distinct image paths, with no conflicting duplicate records. No reference images were downloaded, features regenerated, or replacement sample drawn. Programmatic verification found exactly 942/942 age labels `3-9`; the original has 1,000/1,000. Image-path uniqueness is not proof of unique reference participants.

## Final counts and exclusions

Primary: 1,850 unique known disease patients/images, 379 syndrome labels, 49 main-eligible syndromes and 6,125 syndrome-feature comparisons. Alternative: the same patients and syndrome counts, with 94 different selected photos. Supplementary all-ages: 7,457 patients, 579 labels, 188 main-eligible syndromes.

From 10,942 disease images (8,328 known patients; no missing patient IDs), 9,411 pass frontal quality. Of these, 2,101 meet the recorded childhood age rule; removing 251 additional eligible photos yields 1,850 patients. The other 7,310 frontal-pass images are outside or unresolved for the childhood interval. There are 1,531 quality exclusions. Exclusion reasons are sequential, so a quality-failed image is not also counted as an age exclusion.

Among frontal-pass images: 3,780 zero-zero unresolved; 266 missing/nonnumeric; 11 invalid components; 1,374 under 3; 2,101 ages 3 to under 10; 1,154 ages 10 to under 18; 725 age 18+. The age bands follow the existing automated rules and do not claim independent verification of age semantics.

Broad primary disease groups: African 109, Asian 385, European 939, Others 136, unmapped 281. Detailed groups: African 109, East Asian 114, European 939, Latino/Hispanic 95, Middle Eastern 133, Southeast Asian 14, South Asian 72, unmapped 374. Unmapped patients remain in pooled main analyses but not demographic comparisons. Reference detailed counts: 136, 140, 137, 137, 122, 137, 133 respectively. Full syndrome/group counts, including retained but untested groups, are in `syndrome_group_retention_and_tests.csv`; C endpoint totals count a comparison under both groups and must not be summed as unique tests.

Main retains the ≥10 disease-observation and ≤30% feature-missingness rules. B/C retain ≥2 disease observations per required group and ≥20 reference observations, with existing evidence tiers and priority rules. Small eligible groups are not thereby established as adequately powered.

## Old versus corrected findings

| Analysis / layer | Legacy tests | Corrected tests | Legacy significant | Corrected significant |
|---|---:|---:|---:|---:|
| main / pooled | 27,375 | 6,125 | 8,038 | 811 |
| A_omnibus / layer1_broad_ethnicity | 125 | 125 | 79 | 72 |
| A_pairwise / layer1_broad_ethnicity | 750 | 750 | 239 | 191 |
| B / layer1_broad_ethnicity | 97,500 | 36,000 | 5,929 | 375 |
| C / layer1_broad_ethnicity | 49,500 | 11,500 | 2,623 | 479 |
| A_omnibus / layer2_subcategory | 125 | 125 | 94 | 89 |
| A_pairwise / layer2_subcategory | 2,625 | 2,625 | 817 | 767 |
| B / layer2_subcategory | 97,125 | 33,875 | 5,352 | 452 |
| C / layer2_subcategory | 65,750 | 15,000 | 3,482 | 633 |

Main, B and C significance uses their global BH families; A omnibus uses its feature family and A pairwise uses within-pair BH. These families differ and counts should not be combined. Corrected main significant |g|≥0.8 findings: 379 versus 2,536 in available legacy tables. The manuscript’s 356 syndromes/44,500 tests/9,002 significant findings do not match the available legacy output (219/27,375/8,038). `legacy_reconciliation.csv` retains this discrepancy; the untraced manuscript run is not treated as reproduced.

| Analysis / layer | Legacy significant, now ineligible | Legacy significant, still tested but now nonsignificant | Significant in both | Common tests with wider corrected CI |
|---|---:|---:|---:|---:|
| main / pooled | 5,577 | 1,704 | 757 | 5,796 |
| B / layer1_broad_ethnicity | 1,778 | 3,828 | 323 | 23,764 |
| C / layer1_broad_ethnicity | 1,981 | 495 | 147 | 10,170 |
| B / layer2_subcategory | 1,748 | 3,261 | 343 | 22,237 |
| C / layer2_subcategory | 2,768 | 534 | 180 | 12,803 |

Pairwise A/C comparisons are aligned to lexical group order before joining; signed effects, C intervals, and group-specific fields are reoriented consistently. Original group order is retained for traceability; any original direction-text column still describes its original ordering. This prevents a reordered pair from being mislabeled ineligible. The outer-joined comparison tables explicitly distinguish missing corrected tests from tested nonsignificant findings. For common tests they retain effects, sample sizes, intervals and q-values and report effect and CI-width changes. A larger CI width describes lower numerical precision; legacy versus corrected bootstrap sampling units/methods differ, so it is not a controlled comparison of uncertainty alone. An effect difference or sign reversal is descriptive, not a statistical test of a difference between the two estimates. Loss of significance is not by itself evidence of a changed effect: patient counts, reference deduplication, age restriction, uncertainty estimation and the multiplicity family also changed. All unmatched main hypotheses are excluded by syndrome eligibility in this dataset; B/C tables reflect subgroup eligibility.

Alternative selection produces 834 significant main findings versus 811 primary: 770 are significant in both, 41 primary-only and 64 alternative-only, on the same 6,125 hypotheses. Primary and alternative A use identical references. Full A/B/C sensitivity tables are supplied alongside both main tables; no result was used to reconsider the primary cohort.

## Validation and historical Methods provenance

All 13 unit tests pass. Fresh phenotype-only input exports contain the exact ordered 125 names in `exact_feature_names.csv`, plus only the specified identity/quality/pose fields. The main script’s numeric feature selection selects exactly those names. Demographic masters contain race/reference quality and disease syndrome/ethnicity/subcategory/quality metadata. Patient IDs are nonmissing and unique in selected disease inputs. Main/B/C FDR families, sampled main Mann–Whitney tests, interaction coefficient identities and HC3 p-values were independently validated; A FDR families were also checked. The source scripts themselves still lack comprehensive canonical-name enforcement, so the explicit preflight is a required execution gate.

### Author-confirmed historical FaceKit method

Joshua supplied `e474628f3cb17a83edad59ac6e1d9e5438753032` as the corrected extraction revision after the initially supplied `13603e45add42b7548327373dda673977d2ad08a` could not be resolved. The corrected commit was fetched successfully from the configured HongzhuoChen/FaceKit upstream. Its author and commit date are August 30, 2026. This is author-confirmed extraction provenance with inspected source, not a recovered per-image execution log or model checksum.

The inspected extractor is preserved verbatim as `historical_extractor.py`. Its canonical ordered 125 feature names match the supplied data exactly (`historical_feature_schema.csv`). No image extraction was performed to obtain this schema: only the extractor’s synthetic schema probe was used.

A source-grounded Methods description is:

“Features were extracted with FaceKit at author-confirmed revision e474628f3cb17a83edad59ac6e1d9e5438753032. The extractor operates on two-dimensional landmark coordinates in image pixels; in the image entry point, normalized landmark x and y coordinates are multiplied by image width and height. Landmarks are centered at the midpoint of the inner canthi (indices 133 and 362) and rotated by the negative angle of the line joining them. Distances and areas are normalized within individual feature formulas, principally by bizygomatic width (distance between landmarks 234 and 454) and its square, respectively. Other features use face height, midface height, mouth width, or local ratios; angular features are in degrees. Yaw, pitch and roll are derived from the facial transformation matrix for optional frontal-quality gating. The implementation’s default absolute thresholds are 15°, 15° and 10°, respectively. This is an in-plane roll correction with pose screening, not correction of yaw/pitch foreshortening or three-dimensional frontalization.”

Source locations in this saved revision: landmark definitions lines 40–104; Euler decomposition lines 145–159; threshold defaults lines 185–200; pixel conversion lines 314–320; dropping z in the landmark entry point lines 345–352; optional pose screening lines 377–402; canonicalization/scales lines 411–425; feature formulas lines 457–816. The quality gate is optional and its defaults alone do not prove invocation settings for every historical record. This reanalysis uses the stored `frontal_ok` values and does not reclassify images from default thresholds. Pose-based image ranking uses the existing 15/15/10 denominators, separately from the stored quality decision.

Examples establish why a single universal divisor would be inaccurate: `philtrum_length` uses midface height; `philtrum_width` and lip-height measures use mouth width; `chin_height` uses face height; `nose_to_face_area_ratio` uses a face-area denominator; `cheek_area_r/l` use bizygomatic width squared. Individual functions also implement local ratios, asymmetries and polynomial-derived descriptors.

The local checkout remains `6baac947ac5dbb03ecdbed74615ff53e33460fa6`. After excluding docstrings and formatting, all shared function ASTs in its extractor match the inspected historical extractor, including normalization, Euler decomposition and feature computations. Exact feature names also match. This comparison supports the source description; it does not establish historical detector/model or dependency equivalence. No September 1 implementation was used as a substitute for the author-confirmed revision, and no checkout upgrade or feature regeneration occurred.

## Deliverables and remaining action

`cohorts/` contains local-only row-level selection manifests and corrected exports. `age_3_9/` contains the primary main and A/B/C results; `age_3_9_id_sensitivity/` contains alternative-selection results; `all_ages/` is reused supplementary/exploratory output. Aggregate counts/exclusions, `headline_comparison.csv`, `finding_transition_summary.csv`, `legacy_vs_corrected_*.csv`, and `alternative_selection_comparison.csv` support the conclusions above. Patient-level identifiers remain local; do not distribute the cohort manifests with aggregate reports.

The requested numerical work and historical source inspection are complete. No new age review or primary-cohort decision is required. The package is ready for author review of the results and the proposed Methods wording; the live manuscript and comment replies remain untouched.
