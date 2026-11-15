# Sensitivity output dictionary

All outputs are private aggregate review candidates, including cells with five patients. Exact syndrome
IDs and canonical recorded categories define joins; feature IDs follow the pinned ordered 120-feature schema.
Formulas, anatomical regions and feature-specific units are in exact_120_feature_definitions.csv.

| Output / field | Meaning and units |
| --- | --- |
| analysis_plan.json | Locked hypotheses, feature order/hash, source/code hashes, versions, seeds, thresholds, sizes, nulls and requested repetition goals |
| cohort_audit.json; cohort_candidate_cells.csv | Current source/category/QC exclusions, all syndrome/category cells including zero/ineligible cells, reference capacity, blocked immobility |
| age_sex_exchangeability_audit.csv | Case-versus-other GMDB age absolute SMD (unitless) and maximum recorded-sex fraction gap (0–1); predetermined blocking trigger |
| fairface_placebo_support.csv | Per-cell pseudo-case/reference n, available pool, disjoint design, exact shortage, support flag |
| fairface_placebo_summary.csv | One row per repetition: full-family BH discovery count, by-race counts, minimum supported q, maximum absolute g; primary-count comparison is descriptive |
| gmdb_permutation_aggregate.csv | One row per hypothesis and separate null (within-race, or within-race/sex/age) |
| n_case; n_other; observed_U | Counts of distinct recorded GMDB patients; U counts case/comparator pairwise order with half weight for ties |
| observed_hedges_g | Unitless case-minus-other-GMDB standardized mean difference, exact J correction; undefined for zero pooled variance |
| p_asymptotic_or_exact_auto | Observed MW p using primary auto semantics, not the empirical permutation p |
| q_global_BH_asymptotic_or_exact_auto | Separate descriptive observed-MW family; never substituted for the empirical family |
| null_expected_U | Unrestricted n_case*n_other/2, or exact conditional rank expectation under sex/age block label counts |
| exceedances; B | Number of absolute centered permutation statistics at least observed, including ties; actual number of permutations |
| p_empirical_permutation | (1+exceedances)/(B+1), never zero |
| q_global_BH_empirical | One global planned-family BH adjustment for that empirical null |
| p_value_floor | 1/(B+1); Monte Carlo resolution, not evidence of an exact zero tail |
| q_MC_simultaneous_lower/upper | BH applied to exact binomial tail bounds with 95% simultaneous coverage across that null's Monte Carlo family; conditional Monte Carlo uncertainty only |
| MC_FDR_decision_ambiguous | Lower q<.05 and upper q>=.05; estimated decision remains Monte Carlo-sensitive |
| main_g; main_q | Original syndrome-minus-FairFace effect and primary global q; different comparator, shared cases |
| supported; support_status | Eligibility/comparator support; unsupported planned tests retain p=1 |
| gmdb_permutation_progress.csv | Checkpoint B, null, point discoveries, p floor, zero exceedance count |
| reference_resampling_aggregate.csv | One row per planned hypothesis with original g/q and resampling summaries |
| median_g; g_resampling_q025/q975 | Unitless g median and reference-selection quantiles, not automatically confidence intervals |
| sign_agreement_frequency | Fraction of draws whose g sign agrees with original main g (0–1) |
| significant_frequency | Fraction of draws with full-family BH q<.05 (0–1) |
| median_q; minimum_q; maximum_q | Full-family reference-resampling q summaries; not empirical permutation p-values |
| descriptive_stable | Significant frequency>=.80 and sign agreement>=.95, fixed before outcomes |
| reference_varies | False when a race's reference quota exhausts its pool; g then has no image-selection variability |
| reference_resampling_repetitions.csv | Per-draw total/by-race discoveries, overlap with original discoveries, sign changes and planned family |
| cross_design_robust_candidates.csv | Main discoveries with descriptive reference stability plus significant, direction-concordant empirical results in each enabled GMDB null; not independent replication |
| overall_summary.json | Actual repetition counts, discoveries, overlap, unsupported counts, limitations and completion status |
| checkpoints/ | Private atomic numeric chunks and signatures; reproducible from iteration seeds, never raw patient exports |
| artifact_sha256.json; protected_file_hashes.json | Final artifact and preserved-study hashes; input/execution hashes are in the plan |

Maximum placebo g is descriptive and can be large with sparse pseudo-cases. All q-values control their
separately specified families; they are not clinical probabilities. Patient rows, image paths and random
label assignments are not exported in the review tables.
