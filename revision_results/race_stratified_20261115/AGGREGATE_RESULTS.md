# Experimental redesign results — November 15, 2026

The current real-source 120-feature analysis completed successfully. These are aggregate findings
from the experimental race-stratified study. The validated active study and manuscript remain unchanged.
Sensitivity analyses were deferred at the user's request.

## Cohort and reference

- Current GMDB source: 1,845 selected recorded patients, one quality-passing image per patient,
  conditional age 3 to <10, and 120 finite FaceKit measurements.
- Redesigned eligibility: 736 patients across 67 syndromes and 80 syndrome/category cells,
  each with at least five distinct recorded patients. Counts were rebuilt from current private tables.
- FairFace source: 5,864 accepted age-category 3–9 paths. Excluding four conflicting-label duplicate
  rows and 17 repeated-content copies leaves 5,843 images; 5,049 have a supported canonical category.
- Largest feasible proportional reference under the prespecified 3,409 cap: **889 images**,
  comprising 52 African-, 190 Asian-, and 647 European-annotated images. The European pool limits
  the total. Largest-remainder apportionment, sorted paths, seed `20261115`, and sampling without
  replacement give a maximum absolute composition discrepancy of 0.000477, or 0.0477 percentage points.

The broad crosswalk retains the documented inclusion of FairFace Middle Eastern/West Asian within
Asian. These are harmonized recorded categories, without verified genetic ancestry. Others/unknown
and Latino/Hispanic annotations are excluded from the three-category comparisons. Distinct FairFace
paths or image-content hashes do not establish distinct people.

## Statistical findings

| Analysis | Planned tests | BH q < .05 |
| --- | ---: | ---: |
| FairFace three-category Kruskal–Wallis Analysis A | 120 | 104 |
| Within-category syndrome/reference primary comparisons | 9,600 | 982 |

Analysis A has 16 nonsignificant features. Primary discoveries include 844 with absolute Hedges' g
at least 0.8, with 651 positive and 331 negative effects. All 9,600 planned primary comparisons have
sufficient reference support. Tests are two-sided Mann–Whitney comparisons with one global primary
BH family; g is signed case minus reference and uses pooled sample SD with the exact small-sample
correction. The family differs from the earlier 5,760-test, 955-discovery study, so these counts do
not establish replication or greater accuracy.

## Five-syndrome internal ML benchmark

The count-based selection chooses five different syndromes within the European-annotated group:
Coffin-Siris, Cornelia de Lange, KBG, Sotos and Nicolaides-Baraitser. There are 123 patients overall.
Predictors are the exact 120 FaceKit geometry features. Five-fold outer and three-fold inner
stratified CV separate tuning from held-out prediction, with preprocessing fitted inside training folds.

| Model | OOF balanced accuracy | OOF macro-F1 | OOF accuracy | Top-2 accuracy | Top-3 accuracy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Majority dummy | 20.00% | 0.0846 | 26.83% | — | — |
| Standardized multinomial logistic regression | 39.96% | 0.3918 | 39.84% | 69.92% | 79.67% |
| Standardized RBF SVM | 23.24% | 0.2257 | 26.83% | 57.72% | 69.11% |
| Random forest | 29.62% | 0.2804 | 34.15% | 60.16% | 80.49% |

Balanced chance is 20%. These estimates describe internal nested CV, not independent clinical
validation. SVM and random forest have zero held-out recall for Nicolaides-Baraitser in this run.
Selecting a model from these same scores would require further evaluation. Probabilities are not
clinically calibrated, and confusion counts do not establish shared clinical phenotypes.

## Verification and limitations

All **29 new tests and 47 selected existing regression tests passed**. Verification reconstructed
both BH families, independently checked all 9,600 U statistics, replayed reference sampling and
reconstructed every ML score and confusion matrix from private out-of-fold predictions. Source and
execution-code hashes matched. Protected existing analysis, manuscript and sensitivity files remained unchanged.

Within-category comparisons address recorded-category mismatch. They do not remove the perfect
association of case/reference with image source, hidden demographics, capture artifacts, referral
and selection effects, uncertain age semantics, broad category correspondence, or sparse-cell uncertainty.
FairFace is not clinically certified healthy. These findings do not establish causal anatomy or
diagnostic utility. Sensitivity stability, clinical/HPO comparisons and novelty claims remain unassessed.

The complete report, source/output hash manifests, small-cell tables, feature matrices, patient fold
assignments and predictions remain in ignored local storage. This repository summary contains overall
findings only. [Methods, dependencies and reproduction commands](README.md) accompany the implementation.
