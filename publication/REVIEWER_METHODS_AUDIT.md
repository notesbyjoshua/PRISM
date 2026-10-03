# Reviewer-style methods audit

This review treats the corrected childhood analysis under `revision_results/final_20260914/` as authoritative. It evaluates reproducibility and interpretation, not whether the reported associations are clinically useful.

## Cohort and independence

**RESOLVED - one disease photograph per known patient.** The raw disease table contains 10,942 image records and no missing patient IDs. Among the 9,411 stored frontal-quality passes, 2,101 photographs satisfy the recorded age rule. Those photographs represent 1,850 distinct patient IDs, so 251 additional photographs are removed. The primary selection sorts eligible photographs by normalized absolute pose score and then image ID. The alternative selection changes only the final ordering rule. Both retain one row per patient.

**MANUSCRIPT_LIMITATION - reference independence is not established.** Deduplication removes 58 repeated rows sharing the same file path and leaves 942 paths. A path identifies an image, not a person. The available FairFace metadata do not provide a participant identifier that can establish whether different paths show the same person.

**RESOLVED - childhood age rules are reproducible.** The primary rule is `3 <= age_year + age_month/12 < 10`, with integer year and month components, month 0-11, and explicit exclusion of missing, invalid, and zero-zero values. The same eligibility rule precedes both photograph-selection strategies.

**MANUSCRIPT_LIMITATION - recorded ages are not independently verified.** The code reproduces the supplied values but cannot determine whether every source record uses the same age semantics.

## Multiple testing

**RESOLVED - correction families are distinguishable in the saved outputs.** The pooled main analysis applies one global BH correction across 6,125 syndrome-feature tests. Analysis A omnibus corrections span 125 features per demographic layer; A pairwise corrections are within pair. Analyses B and C each use their own global families, separately by broad and detailed layer. The supplementary tables retain the corresponding q-value names.

**RESOLVED - manuscript examples use the pooled global q column.** The audit selects `mannwhitney_q_fdr_global` for the main examples and `global_fdr_q_value` for matched and interaction analyses. It does not substitute within-syndrome or within-group q-values.

**MANUSCRIPT_LIMITATION - counts from different families are not additive.** Omnibus, matched, and interaction findings answer different questions and use overlapping observations. A total across those tables would have no single multiplicity interpretation.

## Feature dependence

**ADDITIONAL_ANALYSIS_POSSIBLE - dependence is now characterized descriptively.** `publication/feature_correlation/` computes the 125-by-125 Spearman matrix in the 942-path reference cohort and clusters features using `1 - |rho|`. Bilateral measures, mean/asymmetry derivatives, shared landmarks, and common denominators make many tests dependent.

**MANUSCRIPT_LIMITATION - several significant measurements can reflect one facial pattern.** A finding count is not a count of independent biological discoveries. Correlation does not by itself identify identical landmark formulas, and reference-cohort correlation may differ in a syndrome cohort.

**ADDITIONAL_ANALYSIS_POSSIBLE - formula-level dependency map.** The historical extractor could be parsed to list landmarks and denominators used by every feature. That would be more specific than empirical correlation, but it should remain descriptive unless a revised inferential plan is approved prospectively.

## Normalization

**RESOLVED - feature definitions carry normalization context.** Supplementary Table S6 links each exact feature name to its measurement description and formula or normalization context. Philtrum length, for example, is normalized by midface height; several areas use squared bizygomatic width; other measurements use mouth width, face height, or local ratios.

**MANUSCRIPT_LIMITATION - normalized ratios do not isolate the numerator.** A larger ratio can arise from a larger numerator, a smaller denominator, or both. Results should use wording such as “higher normalized philtrum length,” not “anatomically longer philtrum,” unless absolute geometry is analyzed separately.

## Demographic analyses

**RESOLVED - mappings and excluded labels are reproducible.** Existing mapping functions define broad and detailed groups. Unmapped disease patients remain in the pooled analysis but are excluded from demographic comparisons; the saved balance and group-summary tables expose those counts.

**RESOLVED - eligibility thresholds are explicit.** Matched and interaction analyses require at least two disease observations in required groups and at least 20 reference observations, with saved evidence tiers. These are computational eligibility thresholds, not claims of adequate power.

**MANUSCRIPT_LIMITATION - small subgroups are easy to overinterpret.** Some matched effects are estimated from two or three disease observations. Their presence in a complete table is appropriate, but narrative claims should emphasize sample size, interval width, and whether support appears across selections.

**NEEDS_CLINICAL_REVIEW - demographic labels are broad social categories.** Clinical interpretation should avoid treating harmonized race or ethnicity labels as precise genetic ancestry. Review is needed for wording, category aggregation, and any biological explanation.

## Interaction analysis

**RESOLVED - the formal test is the disease-by-group coefficient.** Supplementary Table S4 reports the interaction beta, HC3 standard error, confidence interval, and global q. Group-specific Hedges' g values and delta g are retained as descriptive quantities only.

**MANUSCRIPT_LIMITATION - reversals are unstable when groups are small.** A sign change in descriptive g is not itself an interaction test. Even a significant interaction should be interpreted with the group sizes, interval, multiplicity family, and source limitations visible.

## Photograph-selection sensitivity

**RESOLVED - the comparison is a sensitivity analysis.** The alternative rule uses the same patients, source datasets, eligibility rules, measurements, and reference cohort. It changes 94 selected disease photographs and retains 770 of 811 primary significant comparisons. This is not independent replication.

**RESOLVED - effect changes are exposed.** Supplementary Table S5 includes primary and alternative g, both q-values, direction agreement, and absolute change in g. The overlap figure therefore does not reduce the assessment to significance counts alone.

## Dataset-source confounding

**REQUIRES_EXTERNAL_DATA - disease status is confounded with dataset source.** Disease images come from GMDB and reference images from FairFace. Camera, lighting, compression, cropping, pose, curation, and demographic labeling can therefore differ systematically by disease status. No available metadata field breaks that collinearity, so the analysis cannot estimate a source effect separately from a disease effect.

**ADDITIONAL_ANALYSIS_POSSIBLE - limited sensitivity checks.** Available pose, stored quality, age-bin, sex, and harmonized demographic metadata can support balance descriptions, regression adjustment, or restricted analyses. These checks may show whether results are sensitive to recorded covariates, but they cannot solve unmeasured acquisition differences or perfect source-status confounding.

**MANUSCRIPT_LIMITATION - scope of inference.** The study supports quantitative associations within this assembled cross-dataset comparison. It does not establish diagnostic validity, transportability to clinical photography, or acquisition-invariant anatomy.

## Legacy-result contamination

**RESOLVED - publication generators are pinned to corrected inputs.** The new scripts read `revision_results/final_20260914/` and require 125 exact features, 49 eligible syndromes, and 6,125 main tests. The provenance file records the corrected result and evidence commits.

**MANUSCRIPT_LIMITATION - stale outputs remain in the repository.** They are historical artifacts and must not be used as publication inputs:

- `pds_results_v2/pds_v2_all_disease_feature_results.csv` has 44,500 rows.
- `pds_results_v2/pds_v2_significant_features.csv` has 9,002 rows.
- `pds_results_v2/pds_v2_normalized_matrix_for_ml.csv` and `pds_v2_signed_matrix.csv` each represent 356 syndrome rows.
- `stats_testing/cohort_audit.py:197` hard-codes the legacy tuple containing 579, 356, 44,500, and 9,002 for reconciliation.
- `revision_results/final_20260914/FINAL_REANALYSIS_REPORT.md:39` mentions those values only to document that the manuscript-era run could not be reproduced.

The first four locations can contaminate a new analysis if selected by a broad file glob. The report mention is safe context but may trigger a textual lint warning. Publication scripts should use explicit allow-listed paths, as the new package does.

## Clinical and literature review

**NEEDS_CLINICAL_REVIEW - phenotype wording.** Cornelia de Lange philtrum length and Noonan interpupillary distance are classified here as established clinical associations recovered quantitatively. Rubinstein-Taybi malar bulge remains a candidate underreported quantitative association. Other rows remain unclassified until literature and clinical review are recorded.

**REQUIRES_EXTERNAL_DATA - novelty and validation.** HPO absence cannot establish novelty. Claims of a new phenotype, diagnostic accuracy, or external replication require a documented literature search and an independent clinical or imaging cohort.
