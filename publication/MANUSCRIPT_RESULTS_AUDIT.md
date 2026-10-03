# Manuscript results audit

This audit reconstructs manuscript claims from the corrected historical 125-feature analysis. It does not use newer FaceKit outputs.

**Result:** 90 PASS, 0 FAIL, 0 UNRESOLVED.

| Claim | Expected | Computed | Status | Source and calculation |
|---|---:|---:|:---:|---|
| Disease-image records | 10942 | 10942 | PASS | `data/master_dataset_disease.csv`; row count |
| Stored frontal-quality passes | 9411 | 9411 | PASS | `data/master_dataset_disease.csv`; count frontal_ok values parsed as true |
| Age-eligible frontal-pass photographs | 2101 | 2101 | PASS | `data/master_dataset_disease.csv`; frontal pass and 3 <= age_year + age_month/12 < 10 after validating integer components |
| Repeat-patient photographs removed | 251 | 251 | PASS | `data/master_dataset_disease.csv`; eligible photograph rows minus distinct known patient IDs |
| Retained disease patients | 1850 | 1850 | PASS | `data/master_dataset_disease.csv`; distinct patient IDs among age-eligible frontal-pass rows |
| Retained syndrome labels | 379 | 379 | PASS | `revision_results/final_20260914/syndrome_retention_and_main_tests.csv`; distinct syndrome labels in primary cohort rows |
| Main-analysis eligible syndrome groups | 49 | 49 | PASS | `revision_results/final_20260914/syndrome_retention_and_main_tests.csv`; count main_eligible=true in primary cohort |
| Distinct reference image paths | 942 | 942 | PASS | `data/master_dataset_healthy.csv`; normalized distinct file paths |
| Measurements per eligible syndrome | 125 | 125 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; distinct feature names |
| Primary syndrome-feature tests | 6125 | 6125 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; row count |
| Historical feature schema | 125 | 125 | PASS | `revision_results/final_20260914/exact_feature_names.csv`; row count Notes: Exact schema matches all 125 feature columns used by the result table. |
| Globally significant comparisons | 811 | 811 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count global BH q < 0.05 |
| Significant comparisons with |Hedges' g| >= 0.8 | 379 | 379 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count q < 0.05 and \|g\| >= 0.8 |
| Significant comparisons with |g| < 0.8 | 432 | 432 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count q < 0.05 and \|g\| < 0.8 |
| Nonsignificant comparisons | 5314 | 5314 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count global BH q >= 0.05 |
| Significant overall (%) | 13.2 | 13.240816 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; 100 * 811 / 6125 |
| Significant associations meeting large-effect threshold (%) | 46.7 | 46.732429 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; 100 * 379 / 811 |
| Disease photographs changed | 94 | 94 | PASS | `revision_results/final_20260914/selection_sensitivity_summary.csv`; selected_image_changed |
| Selected disease photographs changed (%) | 5.1 | 5.081081 | PASS | `revision_results/final_20260914/selection_sensitivity_summary.csv`; 100 * changed / 1,850 |
| Alternative-selection significant comparisons | 834 | 834 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count alternative q < 0.05 |
| Shared significant comparisons | 770 | 770 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; set intersection |
| Primary-only significant comparisons | 41 | 41 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; primary set difference |
| Alternative-only significant comparisons | 64 | 64 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; alternative set difference |
| Union of significant comparisons | 875 | 875 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; set union |
| Selection-sensitivity Jaccard overlap (%) | 88 | 88 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; 100 * intersection / union |
| Primary significant findings retained (%) | 94.9 | 94.944513 | PASS | `revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz`; 100 * shared / primary significant |
| Broad reference omnibus significant | 72 | 72 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_A_healthy_kruskal_all_groups.csv.gz`; count corrected omnibus q < 0.05 |
| Broad Analysis B tests | 36000 | 36000 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; row count |
| Broad Analysis B globally significant | 375 | 375 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; count global_significant_fdr=true |
| Broad Analysis C tests | 11500 | 11500 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_C_all_interactions.csv.gz`; row count |
| Broad Analysis C globally significant | 479 | 479 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_C_all_interactions.csv.gz`; count global_significant_fdr=true |
| Detailed reference omnibus significant | 89 | 89 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_A_healthy_kruskal_all_groups.csv.gz`; count corrected omnibus q < 0.05 |
| Detailed Analysis B tests | 33875 | 33875 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_B_matched_disease_vs_healthy.csv.gz`; row count |
| Detailed Analysis B globally significant | 452 | 452 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_B_matched_disease_vs_healthy.csv.gz`; count global_significant_fdr=true |
| Detailed Analysis C tests | 15000 | 15000 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_C_all_interactions.csv.gz`; row count |
| Detailed Analysis C globally significant | 633 | 633 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer2_subcategory/analysis_C_all_interactions.csv.gz`; count global_significant_fdr=true |
| forehead tested comparisons | 245 | 245 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to anatomical region; row count |
| forehead significant comparisons | 74 | 74 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to region; count global q < 0.05 |
| forehead significant (%) | 30.2 | 30.204082 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; 100 * significant / tested |
| chin / jaw tested comparisons | 245 | 245 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to anatomical region; row count |
| chin / jaw significant comparisons | 66 | 66 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to region; count global q < 0.05 |
| chin / jaw significant (%) | 26.9 | 26.938776 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; 100 * significant / tested |
| global face shape / proportions tested comparisons | 294 | 294 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to anatomical region; row count |
| global face shape / proportions significant comparisons | 50 | 50 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; join feature to region; count global q < 0.05 |
| global face shape / proportions significant (%) | 17 | 17.006803 | PASS | `local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results`; 100 * significant / tested |
| Positive significant effects | 490 | 490 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count q < 0.05 and g > 0 |
| Negative significant effects | 321 | 321 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; count q < 0.05 and g < 0 |
| iris_offset_y recurrence count | 24 | 24 | PASS | `corrected main results`; group by exact feature and count q < 0.05 Notes: The manuscript shorthand iris_offset_y corresponds to exact feature iris_offset_y_l. |
| iris_offset_y syndromes tested | 49 | 49 | PASS | `corrected main results`; distinct eligible syndrome IDs for exact feature |
| iris_offset_y recurrence (%) | 49 | 48.979592 | PASS | `corrected main results`; 100 * significant / tested |
| forehead_height recurrence count | 22 | 22 | PASS | `corrected main results`; group by exact feature and count q < 0.05 |
| forehead_height syndromes tested | 49 | 49 | PASS | `corrected main results`; distinct eligible syndrome IDs for exact feature |
| forehead_height recurrence (%) | 44.9 | 44.897959 | PASS | `corrected main results`; 100 * significant / tested |
| hairline_height recurrence count | 22 | 22 | PASS | `corrected main results`; group by exact feature and count q < 0.05 |
| hairline_height syndromes tested | 49 | 49 | PASS | `corrected main results`; distinct eligible syndrome IDs for exact feature |
| hairline_height recurrence (%) | 44.9 | 44.897959 | PASS | `corrected main results`; 100 * significant / tested |
| philtrum_width recurrence count | 21 | 21 | PASS | `corrected main results`; group by exact feature and count q < 0.05 |
| philtrum_width syndromes tested | 49 | 49 | PASS | `corrected main results`; distinct eligible syndrome IDs for exact feature |
| philtrum_width recurrence (%) | 42.9 | 42.857143 | PASS | `corrected main results`; 100 * significant / tested |
| philtrum_length recurrence count | 19 | 19 | PASS | `corrected main results`; group by exact feature and count q < 0.05 |
| philtrum_length syndromes tested | 49 | 49 | PASS | `corrected main results`; distinct eligible syndrome IDs for exact feature |
| philtrum_length recurrence (%) | 38.8 | 38.77551 | PASS | `corrected main results`; 100 * significant / tested |
| Cornelia de Lange syndrome / philtrum_length primary n | 64 | 64 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read n |
| Cornelia de Lange syndrome / philtrum_length primary g | 1.2804 | 1.280422 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Cornelia de Lange syndrome / philtrum_length primary 95% CI lower | 0.9891 | 0.989068 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI lower |
| Cornelia de Lange syndrome / philtrum_length primary 95% CI upper | 1.6064 | 1.606406 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI upper |
| Cornelia de Lange syndrome / philtrum_length primary global q | 7.552e-13 | 7.5521e-13 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read global q |
| Noonan syndrome / inter_pupillary_distance primary n | 34 | 34 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read n |
| Noonan syndrome / inter_pupillary_distance primary g | 1.3189 | 1.318882 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Noonan syndrome / inter_pupillary_distance primary 95% CI lower | 0.82 | 0.820018 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI lower |
| Noonan syndrome / inter_pupillary_distance primary 95% CI upper | 1.8178 | 1.817771 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI upper |
| Noonan syndrome / inter_pupillary_distance primary global q | 2.881e-06 | 2.88068e-06 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read global q |
| Rubinstein-Taybi syndrome / malar_bulge_mean primary n | 24 | 24 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read n |
| Rubinstein-Taybi syndrome / malar_bulge_mean primary g | -1.0016 | -1.001559 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Rubinstein-Taybi syndrome / malar_bulge_mean primary 95% CI lower | -1.4628 | -1.462833 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI lower |
| Rubinstein-Taybi syndrome / malar_bulge_mean primary 95% CI upper | -0.4923 | -0.492294 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read 95% CI upper |
| Rubinstein-Taybi syndrome / malar_bulge_mean primary global q | 0.002809 | 0.002809 | PASS | `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read global q |
| Cornelia de Lange syndrome / philtrum_length alternative g | 1.321 | 1.321004 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Cornelia de Lange syndrome / philtrum_length alternative global q | 2.52e-13 | 2.51958e-13 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read global q |
| Noonan syndrome / inter_pupillary_distance alternative g | 1.32 | 1.316091 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Rubinstein-Taybi syndrome / malar_bulge_mean alternative g | -0.9871 | -0.98706 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read g |
| Rubinstein-Taybi syndrome / malar_bulge_mean alternative CI lower | -1.45 | -1.446622 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read bootstrap CI |
| Rubinstein-Taybi syndrome / malar_bulge_mean alternative CI upper | -0.48 | -0.478074 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read bootstrap CI |
| Rubinstein-Taybi syndrome / malar_bulge_mean alternative global q | 0.003719 | 0.003719 | PASS | `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; select exact syndrome and feature; read global q |
| Cornelia de Lange / philtrum_length / European n | 29 | 29 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read n |
| Cornelia de Lange / philtrum_length / European g | 1.1377 | 1.137662 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read g |
| Cornelia de Lange / philtrum_length / European global q | 0.000834 | 0.000834 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read global q |
| Cornelia de Lange / philtrum_length / Others n | 21 | 21 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read n |
| Cornelia de Lange / philtrum_length / Others g | 1.4387 | 1.438659 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read g |
| Cornelia de Lange / philtrum_length / Others global q | 0.001367 | 0.001367 | PASS | `revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz`; select syndrome, feature, and broad group; read global q |

The script exits nonzero if any checked claim fails or cannot be reconstructed.
