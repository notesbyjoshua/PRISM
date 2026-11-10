# Final targeted 120-feature research completion report

November 9, 2026. Repository base: `2596428d51d423c0e179045a700ab30b43aac1f7`. Active study: `PRISM_120_20261028`.

## Completed scope

Implemented and ran a separate conditional PDIDB–GMDB distributional-consistency analysis on authorized local 120-feature measurements. Added aggregate tables, six original SVG/PNG figures, a source-verified manuscript insertion document, a reassessment of all 17 clinical-literature questions and a FaceKit overlap audit. The local manuscript and its template now distinguish the mentor's generation statement from unverified checkpoint provenance and summarize the separate synthetic results.

The validated primary study remains unchanged: 1,845 selected recorded GMDB patients, 3,409 reference paths, 48 eligible syndromes, 5,760 hypotheses and 955 global-BH associations. No synthetic observations enter primary, demographic, within-GMDB, adjusted or photo/reference inferential families. No extraction or project-wide primary result regeneration was performed. Independent numerical validation recomputed checks against the frozen outputs. No obsolete measurements or comparison results were used in inference.

The Google Doc's Formal Paper tab was inspected read-only to identify placeholders and its 17 appendix questions. `FINAL_120_RESULTS_INSERTIONS.md` is ready for insertion there. No Google Doc tab or comment was changed; BIBM and PRISM were not edited. Public-facing outputs contain only aggregate measurements, category counts, source hashes and original artwork. Private images, individual rows, IDs, encoded-byte hashes and source-patient links were not exported or uploaded.

## Input provenance and roles

| Input / pin | Classification and verified contents |
| --- | --- |
| `ACTIVE_ANALYSIS.json` | Unchanged active-run selector; exact 120 features; zero synthetic primary rows. |
| `<private_run>/COHORTS/primary_master.csv` and `primary.parquet` | Private frozen GMDB selection: 1,845 unique recorded patient IDs and selected image IDs; age 3 <= age < 10; source feature equality verified. Nine matched syndromes contribute 314 selected patients. |
| `<private_run>/COHORTS/reference_balanced_master.csv` and `reference_balanced.parquet` | Private frozen 3,409-path FairFace reference; exact vector equality and source/model pins verified. Same reference used for every unmodified effect profile. |
| `local_review/feature_extraction_120/20261028_current120/gmdb/gmdb_120_all_records_measurements.parquet` | Private current measured source; selected-row geometry, FaceKit/model pins and encoded-byte duplication were checked without displaying records. |
| `<private_run>/PDIDB/pdidb_synthetic_120_{all_records,quality_passed}.parquet` | Private separate synthetic measurements: 1,727 extracted rows; 1,563 QC-passing rows. QC subset equality verified, allowing serialization-only nullable dtype differences. No old public tables substituted. |
| PDIDB revision | `d0907a56525dc3b25712d360f9ce3031150eb535`; source backup SHA-256 `17ec206823ee0c1fefb446173600fb5b9d3ca3c6a84c7fe8983c91529275e980`. Existing extraction provenance read; no credentials accessed and no new dataset download. |
| FaceKit commit | `5764bdd4af9360b55a046a16579fb5a890a8182f`; extractor SHA-256 `8657c0fa5e9aad74d300c8c9bbf023e35f6c5f0a0f74559b2febb43f57fe0a0e`. No source update. |
| MediaPipe model | SHA-256 `64184e229b263107bc2b804c6625db1341ff2bb731874b0bcc2fe6544e0bc9ff`; identical row/source pins across measured datasets. |
| Ordered feature definitions | `publication/feature_definitions/current_120_feature_formulas.csv` and current S6; all 120 ordered semantics matched. Missing/extra/renamed/nonfinite accepted columns fail closed. |

`<private_run>` is `revision_results/rebuild_120_20261028/local_only/run_20261028/`. New private audit reports and input-path/hash classifications are under `local_review/final_research_20261109/`, including `synthetic_consistency/private_input_manifest.json`. They are Git-ignored. Original caches and private read-only archives were not changed. Only archived preparation code was consulted to trace the prespecified adjusted target list; no archived measurement values were loaded for new inference.

The mentor states that StyleGAN was trained on GMDB and PDIDB sampled its learned distribution. The new manifest records this as mentor context received through the November 9 handoff, not provider-verified per-image or checkpoint metadata. The original extraction provenance is preserved. Training membership/counts, model-tag checkpoint identity and equivalence to the FaceKit paper's generator remain unknown.

## Synthetic findings and limitations

The explicit syndrome-label mapping accepts nine exact current names. Generic PDIDB 22q11.2 deletion is excluded because its scope cannot establish equivalence to the current DiGeorge or distal-deletion labels. The matrix also lists all 39 primary cohorts without synthetic counterparts. A blank real n for the unresolved generic label means no equivalent real cohort has been established. No fuzzy matching or source-patient crosswalk was constructed.

The nine supported comparisons use 314 selected GMDB patients and 1,424 synthetic images, with at least 10 rows per source. All five predeclared clinical anchors are available. All 1,727 synthetic records lack age_range, gender and race; equivalent age filtering and demographic adjustment are unavailable. QC-passing tags are v1=734 and v2=829 across the entire synthetic collection. Pooled results describe the observed release mixture, with complete separate tag checks; tags do not establish independent generators.

Across nine equally weighted syndrome summaries, median absolute synthetic-minus-real mean shift is **0.166425 real SD**, median variance ratio **0.986306**, median effect-direction agreement **85.0%**, median Pearson correlation **0.896026** and median Spearman correlation **0.890340**. Each profile uses the identical frozen FairFace reference. Features share geometry and are not 120 independent observations. Fixed-reference-SD contrasts separate mean shifts from pooled-variance changes. No synthetic p-values, FDR family, equivalence conclusion, bootstrap intervals or effective clinical sample-size increase is reported.

Compression is heterogeneous: Cornelia de Lange has median variance ratio 0.838823 and lower synthetic variance in 78.3% of features; KBG's median ratio is 0.784903. Several other syndromes have median ratios above one. The largest exploratory mean discrepancy is KBG/chin_height, +1.247969 real SD, with variance ratio 2.109495. Angelman/inner-canthal spacing and malar geometry also exceed one real SD in absolute mean shift. These are descriptive discrepancies under unmatched age/source distributions, not identified generator mechanisms.

The five anchors retain their effect signs, but agreement in magnitude varies. Real/synthetic Hedges' g is 1.162833/1.208987 for Cornelia de Lange philtrum, 1.243694/1.151991 for Noonan pupil spacing, 1.266137/0.696375 for Williams–Beuren eyebrow area, 1.039430/0.547738 for Williams–Beuren outer-canthal width, and −1.009032/−0.405529 for Rubinstein–Taybi malar geometry.

Model-tag medians differ: v1 shift/variance ratio/Pearson r = 0.192788/0.915234/0.861300; v2 = 0.202423/1.001376/0.893007. Pooled 1%/99% winsorization gives 0.175273/0.983217/0.892617. Low-real-n exclusions and equal-cluster-weighted summaries are supplied. The n >= 20/30 real-cohort cutoffs were added after inspecting initial pooled results and are explicitly exploratory. Available age-restricted GMDB n is only a size proxy; no training-size effect is claimed. No exact duplicate encoded-byte groups or exact synthetic-to-selected-real byte overlaps were found; this is not an identity or privacy assessment.

Methods, source-level audits, all 4,320 syndrome-feature rows across pooled/tag/winsorization strata, cluster summaries and interpretation are in `publication/synthetic_consistency/`. See `synthetic_consistency_results.md`, the predeclared plan and the source/feature/model/output hashes in `synthetic_consistency_manifest.json`.

## Manuscript quantitative completion and scientific audit

- **Recurrence and regions:** `publication/final_research_120/top_five_recurrent_features.csv` reports the top five across all 48 primary cohorts, with denominator, region, direction counts and exact formulas. `regional_discovery_source.csv` is independently regrouped from S1 and asserted against the frozen nine-region table, including 87/240 chin/jaw, 60/192 forehead, 197/960 eyebrows and 35/384 cheeks/malar.
- **Demographic examples:** `demographic_examples.csv` contains 18 support-selected examples: three per A/B/C analysis in each layer, including nonsignificant cells. Selection uses counts and a fixed readable feature sequence, without significance filtering. `demographic_A_group_summaries.csv` preserves the rounded precision encoded in the source group means/medians. Each example has exact source labels, sample counts, formula/units and its correct family-specific q. A omnibus tests are never substituted for B effects or C interactions.
- **Adjustment:** `adjusted_examples.csv` retains all seven prespecified adjusted models and five anchor cells from the full screen, with separate correction families of 4, 3 and 5,880. Reference-adjusted beta values are +0.045197208 for Cornelia de Lange philtrum, +0.035507268 for Noonan pupil spacing and −0.0060628527 for Rubinstein–Taybi malar geometry (exact values in CSV). The latter's age/sex-adjusted within-GMDB result is nonsignificant: beta −0.00220797, HC3 interval [−0.00563492, 0.00121897], four-model q=0.20666045.
- **48 versus 49:** current primary eligibility uses n >= 10. The separate prespecified adjusted list retains FBXW7 with nine current selected patients. The other 48 labels match exactly. Sparse/overlap rules leave FBXW7 and Multiple congenital anomalies-hypotonia-seizures syndrome unsupported in the adjusted screen (120 cells each), while retaining the 5,880 correction denominator. This follows the frozen plan, preparation source and current extension functions; no threshold was relaxed.
- **Comparator differences:** `comparator_overlap.csv` verifies 955 FairFace associations, 222 all-other-GMDB associations, 173 shared, 782 FairFace-only and 49 within-GMDB-only. These are different estimands, not an empirical false-positive rate.
- **Clinical literature:** all 17 appendix questions have unique active S1/S6 matches and current estimates in `clinical_literature_audit_17.csv`. Four support a directly related trait interpretation: Cornelia de Lange philtrum length, Noonan pupil spacing, Wiedemann–Steiner inner-canthal spacing and Sotos forehead width. No source validates the exact normalized ratio. The remaining 13 are retained in the audit as related, scope-mismatched or unverified constructs, not established clinical confirmations. Coffin–Siris long-philtrum support remains unverified in the retrieved source.
- **FaceKit overlap:** `FACEKIT_OVERLAP_AUDIT.md` maps each of the five main anchors and major analysis types to precise preprint sections/figures and current PRISM outputs. The platform, HPO-correspondence question, reference sensitivity and broad synthetic-fidelity question overlap prior FaceKit work. PRISM's specified cohorts, contrasts, correction families and sensitivities constitute its statistical investigation; no new feature extractor or blanket biological novelty is claimed. Exact generator-release equivalence is unresolved.

`manuscript/FINAL_120_RESULTS_INSERTIONS.md` contains titled paste-ready paragraphs, all corresponding tables, figure captions/alt text, a source map and unresolved scientific items. `TABLE_SOURCE_MANIFEST.json` freezes the aggregate input sources for those tables. The main manuscript's earlier validated primary values and design remain intact.

## Figures and visual review

`manuscript/figures/fig5_facekit_prism_workflow.{svg,png}` and `fig6_anatomical_definitions.{svg,png}` are original vector artwork. The anatomical numerator/denominator definitions were checked against S6 and the pinned source, including eyebrow star-polygon area and the signed cheek-to-interpolated-oval malar index. They are schematic, not detector or patient data.

`fig7_regional_discovery_rates.{svg,png}` displays all nine regional significant/tested rates. `fig8_synthetic_anchor_distributions.{svg,png}` shows only Q1/median/Q3 of all five predeclared anchors, with no individual points or extrema. `fig9_synthetic_effect_profiles.{svg,png}` plots syndrome-feature aggregates against an identity line; `fig10_synthetic_fidelity.{svg,png}` shows the nine syndrome summaries. Figures are 7.1-inch layouts, PNG at 300 dpi, and SVG with editable text. Existing Figures 1–4 were unchanged. Captions, alt text and exact source/output hashes are in `FIGURE_CAPTIONS.md` and `FIGURE_SOURCE_MANIFEST.json`.

All six PNGs were inspected. Initial workflow overflow and anatomical arrow/text crowding were corrected and those figures re-rendered and re-inspected. The updated local TeX source compiled successfully using the desktop editor's compiler. This was a source compilation check; no new manuscript PDF was exported or claimed as a deliverable.

## QA commands and evidence

Commands were run from the repository root. The isolated interpreter was `local_review/ci_120_fix/venv/bin/python` (below, `$PY`); it has the pinned public requirements, including PyArrow. The full private validator used the existing `.venv/bin/python`. Analysis commands used `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`; plotting also used `MPLCONFIGDIR=local_review/final_research_20261109/matplotlib`.

```sh
PY=local_review/ci_120_fix/venv/bin/python
$PY stats_testing/synthetic_consistency.py --active-config ACTIVE_ANALYSIS.json --output-dir local_review/final_research_20261109/synthetic_consistency
$PY publication/scripts/final_research_tables.py
$PY publication/scripts/final_research_figures.py
$PY publication/scripts/final_research_text.py
$PY stats_testing/validate_synthetic_consistency.py --report local_review/final_research_20261109/INDEPENDENT_SYNTHETIC_VALIDATION.json
$PY publication/scripts/verify_final_research.py
$PY publication/scripts/verify_reproducible_outputs.py --write-manifest
$PY publication/scripts/verify_reproducible_outputs.py
$PY -m pytest tests -q
$PY -m pytest revision_results/reviewer_followup_2026/test_extended_analysis.py -q
.venv/bin/python revision_results/rebuild_120_20261028/validate.py --scope full
git diff --check
```

Results: **70 repository tests passed** and **16 extension tests passed** (86 total). The sole test warning is an existing Starlette/AnyIO deprecation. New tests cover schema/nonfinite failures, source separation, repeated recorded patients, explicit ambiguous/unmatched mappings, near-zero variance, small n, deterministic computation, fixed-reference signs/scaling, safe exports and source-to-figure/manuscript correspondence. A serialization-only nullable dtype difference between the full and QC PDIDB Parquet files was handled while still asserting equality of values; no observations or measurements were changed.

Independent private synthetic validation checked **108,000 numeric cells across 4,320 rows**, using separate vector calculations and an empirical-CDF integral implementation of Wasserstein distance. Maximum absolute difference was **2.842170943040401e-14**. `INDEPENDENT_SYNTHETIC_VALIDATION.json` records the checked feature-summary hash. Public targeted verification checks all displayed sources, fixed-reference real effects against S1, formulas, support, correction scopes and six figure source/output hashes. The full private validator completed **15 components, status VALIDATED**, with 955/5,760 primary results and zero synthetic primary rows. Frozen config, primary result table, S1 and extractor bytes matched their before-task hashes.

Logs and private evidence: `local_review/final_research_20261109/{synthetic-run,tables,figures,text,independent-synthetic-validation,private-full-validation,tests,extension-tests,package-verification}.log`, `INDEPENDENT_SYNTHETIC_VALIDATION.json` and `primary_before.json`. `publication/ARTIFACT_MANIFEST.json` is intentionally regenerated through the supported verifier after reviewing exports. Public GitHub Actions runs the aggregate checks and synthetic-fixture tests; it does not establish image extraction or patient-level reproducibility. Hosted CI status belongs to the resulting commit's Actions run, separately from these local attestations.

The first hosted run for commit `090d08c` failed at Python setup before installing dependencies or running tests: `ubuntu-latest` had advanced to Ubuntu 26.04, which lacks the pinned Python 3.11.4 binary. The official actions/python-versions manifest lists a Linux x64 3.11.4 binary for Ubuntu 24.04. The workflow is therefore pinned to `ubuntu-24.04`, preserving the Python and scientific dependency pins. This fixes runner selection without changing analyses or loosening checks. See [the setup failure](https://github.com/notesbyjoshua/PRISM/actions/runs/42752136597) and the follow-up commit's run for hosted verification.

## Pending scientific judgments and excluded work

No private-input blocker prevented the requested empirical secondary analysis. Remaining limitations are scientific: clinical confirmation of the nine syndrome mappings; checkpoint identity, training membership and release equivalence; absent synthetic ages/sex/race; the exact clinical meaning of several geometric constructs; and cohort/source mismatch. Broad trait descriptions do not establish novel or clinically validated exact ratios. The analysis uses conditional descriptive point estimates, without claiming uncertainty across generators or new patients.

Google Docs insertion and comment resolution were not performed. No classifier, generator training, diagnostic evaluation or new feature-extraction run was added. Full frontend/browser testing, privacy/IRB/legal publication review and Git-history remediation were explicitly outside this handoff. Existing safeguards and PI-review items remain in force; no release permission was inferred from dataset access or descriptive agreement. The user's unrelated `data/dataset_creation/reloadHFCache.py` remains untouched and excluded from the commit.

## Commit / PR description

**Title:** Add conditional synthetic fidelity and complete current manuscript result support

Add `stats_testing/synthetic_consistency.py` and independent validation, with aggregate PDIDB–GMDB audits/metrics in `publication/synthetic_consistency/`. Add deterministic source-based manuscript tables, literature/FaceKit overlap audits and six original figures through `publication/scripts/final_research_*.py` and `publication/final_research_120/`. Supply paste-ready results and this completion report under `manuscript/`, update the separate synthetic subsection/template, and extend package/source verification and tests. Preserve the frozen 120-feature primary study and keep private inputs excluded. Validation: 86 tests, 108,000 independent synthetic numeric checks, full private validation and aggregate/figure integrity checks.
