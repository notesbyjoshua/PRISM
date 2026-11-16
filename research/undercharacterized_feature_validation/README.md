# Ten locked undercharacterized-feature follow-up candidates

This isolated secondary workstream evaluates internal robustness and follow-up priority, not independent clinical validation or novelty. `analysis_plan.json` fixes thresholds, seeds and the rubric before new patient-level results are calculated. The ten candidate keys and rounded baseline targets are imported from the supplied task markdown, not selected from new outcomes. Exact syndrome/category/feature joins are mandatory.

Set `PRISM_TASK_DOCUMENT` to the authorized task markdown, then run from the repository root:

```sh
.venv/bin/python research/undercharacterized_feature_validation/validate_candidates.py \
  --candidate-document "$PRISM_TASK_DOCUMENT"
.venv/bin/python -m pytest research/undercharacterized_feature_validation/tests -q
```

The default audit input is `local_review/race_stratified_literature_20261115/`; override it with `--audit-directory` only for the same exact locked audit. Missing inputs or baseline mismatches stop the run explicitly. No substitute historical or synthetic data are used. `--output-dir` must pass the existing Git-ignored private-output guard.

The loader follows `SENSITIVITY_ANALYSIS.json` to the experimental race-stratified primary baseline and replays the frozen reference allocation. `ACTIVE_ANALYSIS.json` supplies validated underlying measurements but remains the older manuscript study selector; this work does not promote or edit it. The full 9,600-test primary U/p/g/BH family is independently replayed before candidate analysis, with exact existing tie/method and Hedges correction settings. Source identifiers stay in memory; no photos are accessed.

Outputs stay under ignored `local_review/undercharacterized_feature_validation/`: the ten-row summary in locked order, ranked follow-up list, report, long evidence table, complete cross-category/context tables, private LOO perturbations, aggregate rubric heatmap, frozen plan and hashed manifest. No patient dots or individual distributions are released. The output directory is private, and nothing is copied into public review or committed.

The source tables store normalized features, not raw aligned landmark components. Ratio decomposition is explicitly unavailable; bilateral normalized components are context only. The implementation includes a formula-verification helper for testing component availability and denominator cancellation, but it cannot infer missing numerators. Conditional 2,000-repetition bootstrap intervals resample cases only with fixed references. LOO holds references fixed and reports raw p without a new q. Original q always refers to the full primary family.

Existing 1,000-reference and 100,000-permutation tables are joined with hash checks, not rerun. European reference composition is fixed, so comparator-stability scores are NA there. Within-GMDB permutation effects use a different comparator and conservative simultaneous Monte Carlo decisions remain distinct from point-estimate q. Same-source placebo cannot provide a candidate-specific p-value.

The six-dimensional rubric is descriptive. Completeness and earned/possible scores accompany ranking; NA never becomes zero or an automatic success. Exact formulas and their extraction hash are recorded. Measurement/clinical interpretability scores are conservative heuristics, not clinical grades. Additional targeted review is explicitly required where the prior audit did not include narrative review. No manuscript, cohort, baseline statistics, source datasets or previous sensitivity artifacts are modified.
