# BIBM 2026 evidence report

Generated from the approved childhood PRISM reanalysis. All clinical literature classifications remain pending; no web literature search was performed in this summarization run.

## A. Provenance

- PRISM commit: `a71c9c536a6df42b4c4721668d85ea6234a0be11` (worktree was clean before generation).
- Validated package: `revision_results/final_20260914/`.
- Retained local inputs used only for aggregate cohort and pose checks: `local_review/final_reanalysis_20260914/cohorts/age_3_9/disease_phenotypes.csv` and `cohorts/reference_phenotypes.csv`.
- Historical extraction provenance: `e474628f3cb17a83edad59ac6e1d9e5438753032` and the validated 125-feature schema. No feature extraction was run, and no newer 120-feature FaceKit output or bug-testing result was used.
- Generator: `local_review/bibm2026_evidence/build_bibm_evidence.py`; deterministic, with no new resampling or inferential tests.

Important result-input hashes:

- `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz` — `1e59b26b448e34e5300308cbf15c8cec326efabb1967e8c891d129f9b1672f0f`
- `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz` — `aeaec2e0960a1213169062aa6484f3744e61fe7033ceb90f3903f56d24c4e705`
- `revision_results/final_20260914/exact_feature_names.csv` — `ccbdccc80b350ef3b6e4939e302673bae63931a1958793e9d97d85e76daf8f2a`
- `local_review/final_reanalysis_20260914/cohorts/age_3_9/disease_phenotypes.csv` — `084ee58837e9ea4b76b3b29cca6dbc0313448620be80c4ecb9267c429f11e6f2`
- `local_review/final_reanalysis_20260914/cohorts/reference_phenotypes.csv` — `811e0e5754bbfa0a5dab48e027c754bf97354fdbb4514cbb301f5a7c63994cd5`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_A_healthy_kruskal_all_groups.csv.gz` — `93f13329584eea275d995b2fe1986747c5853c935e075be4275d9c1dacbe889d`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz` — `8757a55b96488b06c111ed0d6f3388934a771b59e99b360a787f02d35a56b1dd`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_C_all_interactions.csv.gz` — `f85759a5630992ccb2e3dcd9721cdcaab1f115db8a59fe2b0d05af4dc83be545`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_A_healthy_kruskal_all_groups.csv.gz` — `d183c54eb2dc760cfcba43fc96c004fcb314ad86efbe019f98fc1e10f2fbb33a`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_B_matched_disease_vs_healthy.csv.gz` — `d3b89732947a2dff5bbc34aa4e6ad6e2a4838a137a5c8401c78b5379de68d251`
- `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_C_all_interactions.csv.gz` — `da7ad0e7073a105cfe1bc1089046d1322991df4920a13bde93040d99fd688c23`

## B. Confirmed headline results

Independent table checks reproduced all frozen counts: 1,850 retained disease patients, 942 distinct reference image paths, 49 eligible syndromes, and 6,125 main tests. Of these, 811 had global BH q<0.05 and 379 also had |g|≥0.8. Primary and alternative photo selections shared 770 significant hypotheses; 41 were primary-only and 64 alternative-only. Analysis A yielded 72/125 omnibus-significant features in the broad layer and 89/125 in the detailed layer.

## C. Paper-ready broader trends

At the exact-feature comparison level, the highest detected-association fractions were forehead: 74/245 (30.2%); chin / jaw: 66/245 (26.9%); global face shape / proportions: 50/294 (17.0%). Across all significant rows, 490 effects were positive and 321 were negative under the stored disease-minus-reference sign convention.

After collapsing right/left/mean variants into transparent feature families and counting each syndrome-family pair once, the leading regional fractions were forehead: 74/245 (30.2%); chin / jaw: 66/245 (26.9%); global face shape / proportions: 50/294 (17.0%). These are detected-association fractions, not formal enrichment tests. Related formulas remain correlated and should not be narrated as independent biological discoveries.

The most recurrent individual feature families were iris_offset_y: 24/49 syndromes (49.0%); forehead_height: 22/49 syndromes (44.9%); hairline_height: 22/49 syndromes (44.9%); philtrum_width: 21/49 syndromes (42.9%); philtrum_length: 19/49 syndromes (38.8%). Recurrence does not imply syndrome specificity or clinical importance.

The complete region, family, syndrome-region, and exact-feature summaries are in the accompanying CSV files. No feature is described as syndrome-specific merely because it was detected once.

## D. Selected syndrome-feature candidates

### 1. Strong quantitative associations needing clinical literature classification

- WILLIAMS-BEUREN SYNDROME; WBS — eb_thickness_mean: g=1.26 (95% CI 1.02 to 1.54), global q=3.15e-13, n=67 [clinical literature classification pending]
- Cornelia de Lange syndrome — philtrum_length: g=1.28 (95% CI 0.99 to 1.61), global q=7.55e-13, n=64 [clinical literature classification pending]
- WILLIAMS-BEUREN SYNDROME; WBS — outer_canthal_distance: g=1.11 (95% CI 0.82 to 1.40), global q=5.89e-11, n=67 [clinical literature classification pending]
- Cornelia de Lange syndrome — philtrum_width: g=-1.01 (95% CI -1.22 to -0.82), global q=8.62e-12, n=64 [clinical literature classification pending]
- Noonan syndrome — inter_pupillary_distance: g=1.32 (95% CI 0.82 to 1.82), global q=2.88e-06, n=34 [clinical literature classification pending]

### 2. Established-looking concordance examples to confirm

These are geometrically interpretable patterns selected for author checking, not claims of established clinical concordance.

- Mucopolysaccharidoses — chin_pointedness_angle: g=1.15 (95% CI 0.86 to 1.44), global q=1.06e-07, n=37 [clinical literature classification pending]
- WIEDEMANN-STEINER SYNDROME; WDSTS — inter_canthal_distance: g=1.88 (95% CI 1.34 to 2.48), global q=1.63e-05, n=14 [clinical literature classification pending]

### 3. Candidates worth investigating as potentially underreported

The label is provisional and does not imply novelty.

- MANNOSIDOSIS, ALPHA B, LYSOSOMAL; MANSA — inter_pupillary_distance: g=1.54 (95% CI 1.02 to 2.13), global q=0.000244, n=14 [clinical literature classification pending]
- Rubinstein-Taybi syndrome — malar_bulge_mean: g=-1.00 (95% CI -1.46 to -0.49), global q=0.00281, n=24 [clinical literature classification pending]

## E. Demographic and interaction examples

Broad and detailed demographic layers are reported separately and are not independent replications. Small nonsignificant subgroups are not treated as evidence of no effect.

- Cornelia de Lange syndrome — philtrum_width, broad/European: n=29 disease and 137 reference; g=-1.34 (95% CI -1.67 to -1.06); global q=5.72e-06; primary; pooled direction agrees.
- Cornelia de Lange syndrome — philtrum_length, broad/European: n=29 disease and 137 reference; g=1.14 (95% CI 0.65 to 1.62); global q=0.000834; primary; pooled direction agrees.
- WILLIAMS-BEUREN SYNDROME; WBS — eb_thickness_mean, broad/Asian: n=25 disease and 532 reference; g=1.32 (95% CI 0.87 to 1.77); global q=0.000295; primary; pooled direction agrees.

Formal high-confidence interaction examples matching the shortlist:

- Cornelia de Lange syndrome — philtrum_width, detailed Latino/Hispanic versus European: interaction beta=-0.0114 (HC3 95% CI -0.0182 to -0.00447), global q=0.0334; descriptive g values -0.52 and -1.34.

## F. Photograph-selection sensitivity

All 25 shortlisted effects had the same direction, were globally significant in both selections, and had confidence intervals excluding zero in both selections by construction. The median absolute change in g was 0.003; the maximum was 0.121. This is a photograph-selection robustness analysis, not independent replication.

## G. Pose-flag diagnostic

- retained disease: 1,850 within thresholds; 0 outside at least one threshold; 0 missing/nonfinite pose.
- retained reference: 942 within thresholds; 0 outside at least one threshold; 0 missing/nonfinite pose.

The comparison uses historical defaults (|yaw|≤15°, |pitch|≤15°, |roll|≤10°) only as a consistency reference. It does not establish the invocation settings used for every historical record, does not override stored `frontal_ok`, and did not change the cohort.

## H. Figure recommendations

1. Use `figure_cross_syndrome_heatmap.png` as the primary Results figure. It uses the ten largest eligible syndrome groups with n≥20 and twelve representative exact features selected from recurrent/large effects while taking one interpretable representative per family and spanning all anatomical regions. Significant cells show signed Hedges' g; tested nonsignificant cells are gray; untested cells would be black (none occur in this complete selected grid).
2. If space permits, derive a small forest panel from the top 6–10 rows of `main_candidate_shortlist.csv`; the stored confidence intervals and alternate-selection columns support this without new inference.

Proposed heatmap caption: **Cross-syndrome structure of selected facial measurements in children ages 3–9.** Columns are the ten largest eligible syndrome groups (n≥20); rows are twelve representative exact measurements from distinct feature families, selected using recurrence, large-effect frequency, interpretability, and anatomical coverage rather than q-value alone. Color encodes signed Hedges' g only for global-BH-significant comparisons (q<0.05); gray cells were tested but nonsignificant. Positive values indicate larger measurements in the syndrome group than in the reference set. Feature variants are not independent, and the display is descriptive.

## I. Exact manuscript replacement text

### Selected phenotype associations and anatomical trends

Among 6,125 syndrome-feature comparisons, 811 met the global Benjamini–Hochberg threshold (q<0.05), including 379 with |Hedges' g|≥0.8. The largest detected-association fractions at the exact-comparison level were forehead: 74/245 (30.2%); chin / jaw: 66/245 (26.9%); global face shape / proportions: 50/294 (17.0%). Family-level de-duplication retained the same broad emphasis (forehead: 74/245 (30.2%); chin / jaw: 66/245 (26.9%); global face shape / proportions: 50/294 (17.0%)), although correlated measurements should not be interpreted as independent biological findings. Representative stable associations included WILLIAMS-BEUREN SYNDROME; WBS eb thickness mean (g=1.26, 95% CI 1.02 to 1.54, global q=3.15e-13, n=67); Cornelia de Lange syndrome philtrum length (g=1.28, 95% CI 0.99 to 1.61, global q=7.55e-13, n=64); WILLIAMS-BEUREN SYNDROME; WBS outer canthal distance (g=1.11, 95% CI 0.82 to 1.40, global q=5.89e-11, n=67); each remains **[clinical literature classification pending]**. All 25 shortlisted effects retained direction, global significance, and confidence intervals excluding zero under alternative eligible-photograph selection; the median absolute change in g was 0.003. This sensitivity analysis supports robustness to the prespecified photograph-selection rule but is not an independent replication.

## Quality-control record

- Main candidates were filtered with `mannwhitney_q_fdr_global`; within-syndrome q-values were not used for selection.
- Primary and alternative rows were joined one-to-one on the exact stored disease label and feature, then labels were cleaned only for presentation.
- All stored main and alternative confidence intervals had correctly ordered bounds; sign agreement was computed directly from signed Hedges' g.
- Denominators use actual tested rows. Untested status is distinct from nonsignificance; the selected main grid happens to be complete.
- Analysis B layers remain separate. Analysis C output uses the formal interaction beta and HC3 interval/q; delta-g is retained only as descriptive context.
- The output directory contains no patient IDs, patient-level manifests, or images. The pose table contains aggregate counts only.
- All outputs derive from the saved historical 125-feature result schema; no updated FaceKit code or regenerated feature values were used.
