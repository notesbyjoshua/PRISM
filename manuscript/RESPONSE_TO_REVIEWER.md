# Response to reviewer

This response is based on the corrected analysis in `revision_results/final_20260914/` and the historical 125-feature FaceKit implementation at `e474628f3cb17a83edad59ac6e1d9e5438753032`. No newer FaceKit result and no new image-level feature extraction is used. Statistical claims below are limited to the verified manuscript-results and methods audits.

## 1. Repeated patients and images

### Reviewer concern

Multiple photographs from the same disease patient could be treated as independent observations and inflate the effective sample size.

### What was changed

The childhood disease cohort now retains one eligible photograph per known patient identifier. After the quality and age filters, 2,101 eligible photographs represented 1,850 patients; 251 additional photographs were removed by the one-image-per-patient rule. The same patient set and eligibility rules are used for the primary and alternative photograph selections.

### Analysis/evidence

The primary selection ranks a patient's eligible photographs by the maximum normalized absolute pose score, `max(|yaw|/15, |pitch|/15, |roll|/10)`, and then lexical image ID. The alternative selects by lexical image ID after the same eligibility filters. The alternative changed 94 selected photographs. Of 811 primary global-BH-significant comparisons, 770 remained significant under the alternative selection.

### Manuscript change

The Cohort construction subsection now states the one-image-per-patient rule, the exact primary ranking, the alternative rule, and why the comparison is a sensitivity analysis rather than independent replication. Results report the eligible-image and retained-patient counts and the overlap of significant hypotheses.

### Remaining limitation

Photograph selection can reduce repeated observations but cannot show that either selected image is a patient's uniquely representative appearance. The sensitivity analysis reuses the same patients and is not external validation.

## 2. Age mismatch

### Reviewer concern

Disease and reference ages may not be comparable, and ambiguous age records could enter the childhood analysis.

### What was changed

Disease eligibility is now stated exactly as `3 <= age_year + age_month/12 < 10`, with integer year and month components, month 0--11, and exclusion of missing, invalid, and zero-year/zero-month records. FairFace references were restricted to the supplied `3-9` category.

### Analysis/evidence

Among disease photographs with a passing stored quality flag, 2,101 met the recorded childhood rule. All 942 retained reference paths carried the FairFace `3-9` label. The audit reproduces these counts from the frozen inputs.

### Manuscript change

The Methods distinguish the GMDB year/month calculation from the FairFace categorical label and state that recorded ages and cross-dataset age semantics were not independently verified.

### Remaining limitation

Exact age matching is impossible because the reference data provide an age category rather than exact age. The analysis also cannot establish that source records use identical age semantics.

## 3. FairFace sampling and deduplication

### Reviewer concern

Repeated FairFace rows or photographs could inflate the reference sample, and distinct files may not represent distinct people.

### What was changed

Reference rows were deduplicated by normalized image path. Fifty-eight repeated-path rows were removed from 1,000 rows, leaving 942 distinct image paths; no replacement reference sample was drawn.

### Analysis/evidence

The reviewer-methods audit and manuscript-results audit independently reproduce the 942-path count. They also distinguish file-path uniqueness from participant uniqueness.

### Manuscript change

The Cohort construction subsection now describes path deduplication and states that the available FairFace metadata do not contain a participant identifier capable of establishing unique people.

### Remaining limitation

Different image paths may show the same person. Reference-observation independence therefore remains unproven.

## 4. Syndrome eligibility

### Reviewer concern

The manuscript did not make clear which syndrome groups were tested or how small groups were excluded.

### What was changed

The primary analysis requires at least 10 retained disease patients per syndrome and retains the existing feature-level requirement of no more than 30% missingness. Forty-nine syndromes met the main-analysis criterion, yielding 49 x 125 = 6,125 primary syndrome--feature tests. The subgroup analyses retain their separate minimums of two disease observations in required groups and 20 reference observations.

### Analysis/evidence

The frozen syndrome-retention table identifies 49 main-eligible groups, and the main result table contains exactly 125 features for each eligible syndrome. Saved subgroup tables expose their cell counts and evidence tiers.

### Manuscript change

The cohort and statistical Methods now give the primary and subgroup eligibility thresholds and distinguish computational eligibility from adequate power.

### Remaining limitation

Some eligible demographic subgroup estimates use only two or three disease observations. Their inclusion in complete tables does not make them stable or adequately powered.

## 5. Feature normalization

### Reviewer concern

The normalization and geometric meaning of the 125 facial measurements were not reproducible from the manuscript, and ratio changes risked being described as absolute anatomy.

### What was changed

The Methods are pinned to historical FaceKit revision `e474628f3cb17a83edad59ac6e1d9e5438753032`. Supplementary Table S6 records every historical expression, source location, landmark dependency, unit, and denominator. The manuscript now describes pixel conversion, inner-canthal centering and roll alignment, and the feature-specific normalization scheme.

### Analysis/evidence

All 125 corrected-result feature names map to one exact source-derived formula. Examples include philtrum length divided by midface height, chin height divided by face height, and cheek area divided by squared bizygomatic width. No newer 120-feature output is substituted.

### Manuscript change

The Historical feature extraction subsection now explains the exact provenance and warns that a normalized ratio can change through its numerator, denominator, or both.

### Remaining limitation

The saved measurements do not support absolute-anatomy claims, three-dimensional frontalization, or reconstruction of a per-image historical execution log. [HUMAN_INPUT_REQUIRED: the FaceKit author should review the exact S6 before submission.]

## 6. Broad and detailed demographic layers

### Reviewer concern

The relationship between broad and detailed demographic analyses was unclear, creating a risk that agreement would be described as replication.

### What was changed

The broad layer is identified as a harmonization that provides larger cells and cross-source comparability; the detailed layer preserves more source-label granularity. Both layers are reported because they expose the sensitivity of results to grouping resolution.

### Analysis/evidence

Analysis A found 72 of 125 omnibus-significant features in the broad layer and 89 of 125 in the detailed layer. Analyses B and C use separate global BH families for each layer. The layers reuse observations and related mappings.

### Manuscript change

The Statistical analysis subsection explains the purpose of each layer, their separate correction families, the exclusion of unmapped records from demographic comparisons, and why the layers are not independent replications.

### Remaining limitation

The harmonized labels are broad dataset annotations, not precise genetic ancestry. Small cells and source-specific labeling can make results unstable or difficult to interpret. [HUMAN_INPUT_REQUIRED: clinical and demographic-label wording requires author review.]

## 7. Interaction interpretation

### Reviewer concern

Differences or sign reversals between subgroup-specific Hedges' g values were at risk of being interpreted as formal interaction tests.

### What was changed

The formal test is now explicitly the disease-by-group coefficient from Analysis C, with HC3 robust standard error, confidence interval, and the Analysis C global q-value. Group-specific Hedges' g and delta-g are retained only as descriptive quantities.

### Analysis/evidence

The saved interaction tables contain the coefficient, HC3 uncertainty, and separate broad- and detailed-layer correction families. Independent validation reproduced the interaction coefficient identities and HC3 p-values.

### Manuscript change

The Methods state that delta-g and sign reversal are not substitutes for the interaction coefficient. Results report Analysis C totals by layer without adding them across correction families.

### Remaining limitation

Even a significant interaction can be unstable when a disease subgroup is small. Interpretation must retain the group sizes, interval width, multiplicity family, and cross-dataset source limitations.

## 8. Candidate phenotype interpretation

### Reviewer concern

Large or recurrent quantitative associations could be presented as novel or clinically established phenotypes without a verified literature review or clinical assessment.

### What was changed

The expanded Results report only audited numerical findings and do not assign clinical meaning to selected associations. Candidate evidence remains a review queue rather than manuscript validation, and HPO absence is not used as evidence of novelty.

### Analysis/evidence

The candidate table records quantitative effects, sensitivity results, and unresolved review status. The BIBM evidence report states that clinical literature classifications remain pending. The manuscript-results audit verifies the numerical headline values but does not verify clinical interpretation.

### Manuscript change

The Results captions and narrative are descriptive. The Discussion and References retain explicit markers rather than filling unverified clinical or literature claims.

### Remaining limitation

[HUMAN_INPUT_REQUIRED: complete a documented literature search and clinical review before classifying associations as established, underreported, or novel, or before making diagnostic claims.]

## 9. Source-dataset confounding

### Reviewer concern

Disease--reference differences may reflect GMDB-versus-FairFace acquisition and curation rather than disease-associated anatomy.

### What was changed

The manuscript now identifies GMDB and FairFace directly and states that disease status is perfectly confounded with dataset source. Three predefined recorded-covariate sensitivity models adjust highlighted comparisons for sex label, broad demographic group, yaw, pitch, and roll, with HC3 uncertainty and BH correction across the three models.

### Analysis/evidence

The common metadata support balance descriptions and limited adjustment, but no field separates a disease effect from a source effect. Exact reference age is unavailable, and unmeasured camera, lighting, compression, cropping, curation, and labeling differences remain.

### Manuscript change

Dataset-source confounding is now stated in Study design, Statistical analysis, and Limitations. The covariate models are described as exploratory sensitivity checks rather than replacements for the primary analysis.

### Remaining limitation

Recorded-covariate adjustment cannot solve perfect source-status confounding or unmeasured acquisition differences. External validation in a source-balanced or independently collected cohort remains necessary for acquisition-invariant or clinical claims.

## Items requiring human input before submission

- **HUMAN_INPUT_REQUIRED:** verified literature citations and clinical-review decisions.
- **HUMAN_INPUT_REQUIRED:** institution-approved ethics, consent, exemption, and oversight language.
- **HUMAN_INPUT_REQUIRED:** funding, conflicts-of-interest, and author-contribution statements.
- **HUMAN_INPUT_REQUIRED:** mentor-, institution-, and data-provider-approved data-release permissions and availability language. The repository is currently private; public redistribution remains pending those decisions.
