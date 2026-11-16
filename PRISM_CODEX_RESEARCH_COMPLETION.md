# PRISM high-value research completion — November 10, 2026

Completed the available, authorized computational work. The most defensible main narrative is **comparator and normalization sensitivity in FaceKit-derived facial geometry associations**. This is exploratory observational research, without a new extractor, clinical validation or validated ranking algorithm.

## Verified baseline and preservation

Baseline Git HEAD: eef1cc0cd4c2be780abdfeae2fb185c4529ee130. Sole active study: PRISM_120_20261028.

- FaceKit commit: 5764bdd4af9360b55a046a16579fb5a890a8182f.
- Extractor SHA-256: 8657c0fa5e9aad74d300c8c9bbf023e35f6c5f0a0f74559b2febb43f57fe0a0e.
- Model SHA-256: 64184e229b263107bc2b804c6625db1341ff2bb731874b0bcc2fe6544e0bc9ff.
- 120 ordered features; 1,845 selected unique recorded GMDB patients; 872 in 48 primary cohorts; 3,409 FairFace reference paths.
- 5,760 primary tests; 955 FairFace-significant, 222 other-GMDB-significant, 173 shared, 782 FairFace-only, 49 GMDB-only, 4,756 neither.
- Zero PDIDB observations in primary or new real-patient inference.

The 344 pre-existing tracked/user-owned files and protected private input sources were hash-checked throughout. Frozen results, ACTIVE_ANALYSIS.json, FaceKit/model pins, public manifest, dashboard, manuscripts and the user's reloadHFCache.py remain unchanged. No historical measurements were substituted. No image pixels were processed, patient data uploaded, tokens accessed, commits/pushes made or history rewritten.

## Deliverables

New code: stats_testing/robustness_20261110/. Tests: tests/test_robustness_20261110.py. Exact commands are in the module README.md. All new empirical outputs remain **private** in:

local_review/prism_robustness_20261110/

| Milestone | Exact private deliverables |
| --- | --- |
| A: preflight/comparators | 00_preflight.md, 00_baseline_validation.json, 00_work_plan.md, MILESTONE_A_CHECKPOINT.md, evidence/P0_exhaustive_evidence_matrix.csv, P0_comparator_evidence.md, P0_source_bias_audit.md |
| B: source/dependence/stability | P1_feature_dependence.md, P1_stability_report.md, P1_PHOTO_REPLAY_VALIDATION.json; controls, proxy distributions, restrictions/models, clusters/edges and photo/patient/estimator/influence summaries in aggregate_candidates/ |
| C: geometry/candidates | P1_clinical_geometry_case_studies.md, P1_ranked_candidate_evidence.csv, P1_selection_rules.md; alternative ratios, numerical simulation, seed/age/family/cluster sensitivities |
| D: scientific decision/literature | P2_method_decision.md = NOT_JUSTIFIED; unconditional literature work completed in P2_phenotype_source_matrix.csv, P2_literature_claim_audit.md, P2_clinical_comparison_manuscript_draft.md, PUBLIC_SOURCE_NOTES.json |
| E: author review | P3_editorial_decision_memo.md, P3_manuscript_insertions.md, P3_figure_plan.md, P3_claim_to_evidence_ledger.csv, P3_reporting_checklist.md; 15 figure drafts in figures/ with editable SVG/PNG and source/code/output hashes |

Every primary hypothesis is joined by canonical syndrome ID and exact feature ID. Separate long tables retain targeted 22,320-cell and age 11,760-cell families; the adjusted 5,880 screen and smaller prespecified families remain distinct. Unsupported estimates/display q remain missing, not proven nulls.

## Findings and scientific assessment

- **Comparator dependence:** 4,107/5,760 effects (71.3%) have the same direction. Among 482 primary q<0.05, |g|≥0.8 associations, 211 (43.8%) more than halve against other GMDB conditions; nine reverse sign. The designs share cases and answer different questions, so overlap is not independent replication or a clinical false-positive estimate. Source: P0_comparator_evidence.md and exhaustive matrix.
- **Acquisition susceptibility:** 100 stratified random FairFace control repetitions yield no discoveries; deliberate within-stratum pitch imbalance yields median 25, range 0–73, correcting 5,760 tests per repetition. This demonstrates workflow sensitivity, without separating disease from source. Source: aggregate_candidates/P0_negative_control_repetitions.csv.
- **Redundancy:** each global source yields 55 descriptive correlation groups; adjusted Rand agreement is 0.857. Groups are not independent biological dimensions; the original testing denominator is unchanged. Source: aggregate_candidates/P1_cluster_source_comparison.csv.
- **Limited photo opportunity:** only 54/872 primary photos change under the deterministic alternate rule; 94/872 patients have multiple eligible images. Exact recovery of 100 existing random-photo runs matches seven frozen summaries across 5,760 pairs, maximum absolute difference 3.2863e-14. Source: P1_PHOTO_REPLAY_VALIDATION.json.
- **Definition sensitivity:** Williams–Beuren outer-canthal span against FairFace changes from g=+1.039 under bizygomatic normalization to −0.103 under height normalization. Cornelia de Lange upper-lip-peak span changes from −0.399 to −0.342. Denominators change the estimand; these do not establish absolute anatomy, a robust reversal or a superior feature. Source: aggregate_candidates/P1_normalization_cases.csv, a new separate 12-test BH family.
- **Positive and weakening examples:** Cornelia de Lange philtrum length, Noonan pupil spacing and Williams–Beuren eyebrow area remain positive against both comparators. Rubinstein–Taybi malar geometry attenuates from primary g=−1.009 to other-GMDB g=−0.245 (other-GMDB global q=0.561). Eyebrow area is not hair density; lip peaks need not be philtral columns. Exact outputs/correction scopes: P3_manuscript_insertions.md.
- **Covariate/age sensitivity:** Williams–Beuren eyebrow area and outer span remain positive but do not meet q<0.05 in the new 12-model quality-adjusted family. All six age-case signs agree across strata, with no younger-stratum and three older-stratum discoveries under the original 11,760-cell age correction. Different support/composition prevents interpreting this as a formal age interaction. Source: P1_covariate_age_sensitivity.md.
- **Uncertainty:** new case-only fixed-reference nominal-95% bootstrap ranges cover 89.0%–94.3% across four numerical regimes, 300 independent datasets each. These remain descriptive and do not automatically invalidate frozen two-sample intervals. No biased concordance interval is reinstated. Source: aggregate_candidates/P1_bootstrap_coverage_simulation.csv.
- **Exploratory selection:** magnitude top-20 lists retain 16–19 pairs across same-patient photo draws. Three extra seeds retain all 20 original high-stability candidates, but this is Monte Carlo sensitivity on the same patients, not held-out performance. Family/cluster collapsing changes other rankings. Source: ranking sensitivity tables and P1_selection_rules.md.
- **Clinical novelty unproven:** all 17 prior numerical literature questions match current outputs, with no exact normalized clinical-measure equivalence established. The 120 potential construct mappings are not 120 validated HPO equivalences. Other candidate mappings are explicitly unverified. Full-text/indexed/indirect/abstract-only access states are distinguished; absent terminology is never evidence of novelty.

## Validation actually performed

Used local_review/ci_120_fix/venv/bin/python, PYTHONPATH=.:FaceKit/src and one-thread BLAS/OpenMP settings. Full exact commands and log locations are in the code README and private validation/ directory.

| Command/check | Actual result |
| --- | --- |
| python -m pytest tests -q --tb=short | **79 passed**, one existing dependency deprecation warning |
| python -m pytest revision_results/reviewer_followup_2026/test_extended_analysis.py -q --tb=short | **16 passed** |
| python publication/scripts/verify_reproducible_outputs.py --report local_review/prism_robustness_20261110/validation/public_aggregate_verification.md | **PASS**, unchanged committed aggregate package; no manifest-writing flag |
| python -m stats_testing.robustness_20261110.validate --full-private-baseline | **PASS**, unchanged frozen validator replayed all 15 private components; only two final report writes redirected into new private output |
| Final python -m stats_testing.robustness_20261110.validate | **PASS**, 11,520 point estimates, 100 photo BH families, new geometry/quality BH, two control replays, aggregate schemas, 15 figure records, 14 claim sources and three seed summaries |
| Visual review | All 15 drafts inspected; crowded labels and an acquisition-plot source-name mismatch corrected before finalization |
| python -m stats_testing.robustness_20261110.provenance --check | Final source/output/input/baseline hash audit; private manifest RESEARCH_PROVENANCE.json |

The first control allocation incorrectly retained pitch order between pseudo-cohorts. Its invalid v1 checkpoint is preserved privately for audit; final results exclusively use corrected v2 allocation with a regression fixture. An exact-token dependency-parser correction also has a fixture and fresh v2 cluster caches. Neither changed primary results.

These are actual local numerical and aggregate checks. No new code was pushed or remote CI run claimed. Earlier CI success applies only to the committed baseline.

## Blocked or unsupported work

1. **BLOCKED_PRIVATE_INPUT:** compatible saved landmarks and separate image-space numerator/denominator components were unavailable. Normalized ratios cannot recover them uniquely.
2. **BLOCKED_IMAGE_REPROCESSING_APPROVAL:** real-image crop/brightness/resolution/pose perturbations were not performed. Compatible inputs and confirmed permission are needed.
3. **NOT_JUSTIFIED:** no objectively evaluated advantage supports a new score. Conditional method specification/benchmark/ablation deliverables were intentionally not implemented. The decision memo defines possible future targets, baselines, leakage safeguards and falsification criteria.
4. FairFace exact ages/health/unique identities, native resolution/crop/confidence/exposure and an independent source-matched clinical reference remain unavailable. Adjustment cannot resolve source/disease confounding. PDIDB does not resolve it.
5. PI review must settle protocol/consent, image-processing/publication rights, rare-cohort disclosure, attribution, historical remote exposure, clinical wording, authorship/funding and release statements. No approvals were invented.

## Prioritized next actions

Privately review the editorial memo, claims ledger and figure plan with the author/PI. Preserve the comparator/measurement sensitivity narrative and unfavorable findings. Clinical review should assess exact measurement comparability. Only pursue component/image experiments when their inputs and permission are confirmed.

Manuscript incorporation, public tables/figures, public-manifest updates, submission and pushing require the separate authorization specified in this handout. No active manuscript or live Google Doc was overwritten.
