# Race-stratified head-pitch imbalance sensitivity

This isolated experiment compares random and deliberately pose-selected pseudo-cohorts drawn only from the accepted, deduplicated FairFace reference pool of `PRISM_RACE_STRATIFIED_20261115`. GMDB supplies the locked eligible cell counts; its facial measurements never enter either pseudo arm. The primary baseline, selectors, previous sensitivities, released figures, manuscript and publication materials are read-only.

## Run

Use the existing project `.venv` from the repository root; no image processing, network access or package installation is required.

```sh
.venv/bin/python research/race_stratified_pitch_sensitivity/run_pitch_imbalance.py --dry-run
.venv/bin/python -m pytest research/race_stratified_pitch_sensitivity/tests -q \
  --junitxml=local_review/race_stratified_pitch_sensitivity/tests.xml
.venv/bin/python research/race_stratified_pitch_sensitivity/run_pitch_imbalance.py \
  --repetitions 10 --output-dir local_review/race_stratified_pitch_sensitivity_smoke
.venv/bin/python research/race_stratified_pitch_sensitivity/validate.py \
  --output-dir local_review/race_stratified_pitch_sensitivity_smoke
.venv/bin/python research/race_stratified_pitch_sensitivity/run_pitch_imbalance.py --repetitions 1000
.venv/bin/python research/race_stratified_pitch_sensitivity/validate.py
```

The CLI also accepts `--seed`, `--output-dir` and `--resume`. Defaults are seed `2026111704`, 1,000 repetitions per arm and ignored `local_review/race_stratified_pitch_sensitivity/`. Dry-run validates sources and prints support without running outcomes or creating a plan. Missing/invalid sources produce a specific BLOCKED audit. Runs below 1,000 are explicitly PILOT, including after verification. A full run remains awaiting independent validation until `validate.py` succeeds.

Plans are frozen before outcomes. Existing plans require explicit `--resume`; source, support, numerical code, library versions, seed and repetition budget must match. Changed plans or numerical code require a fresh ignored directory. Atomic checkpoints cover 25 paired repetitions (or the remaining final block), record SHA256 hashes and enforce contiguous ranges. Resume reuses completed chunks without changing random draws. Reporting code is hashed separately; same-plan resume can regenerate reports without numerical recomputation. Never resume or replace a study using outcomes to choose a design.

## Locked sources and pose provenance

`design.py` reuses the existing validated `load_baseline()` loader, including primary manifest/artifact/source/code hashes, reference selection replay, age/QC rules, one-image-per-patient selection, content deduplication and the broad recorded-category crosswalk. It checks the 67-syndrome, 80-cell, 736-patient, 120-feature, 9,600-test baseline and 982 original BH findings. These baseline values are validation requirements, not results of the new placebo experiment.

Pose columns must be finite and agree exactly with the accepted original FairFace extraction, linked through the verified `image_path_key`/`file` key and content hash. The extractor bytes must match the extraction provenance and pinned FaceKit commit. Pitch is the signed small-angle XYZ Tait–Bryan estimate from the MediaPipe transformation matrix, converted to degrees; the regular branch is `atan2(R[2,1], R[2,2])`. Absolute pitch measures deviation from the extractor zero. Physical up/down interpretation and individual pose error are not independently verified. Existing accepted pose limits are read from extraction provenance, never tuned using discovery outcomes. Nonfinite, missing, unverified or QC-inconsistent pose blocks execution; exclusions are not silently substituted.

The availability audit includes signed/absolute pitch, yaw and roll distributions, missingness and ties. The current extractor already pose-canonicalizes features; this experiment tests selection sensitivity in those existing normalized measurements without changing their definitions.

## Prespecified design

Within each category, sort absolute pitch ascending, breaking ties by the verified lexical private file key. The lower `floor(N/2)` rows form the low half; the odd row belongs to the high half. Each cell must support its exact original case count and reference quota in both arms. Freeze their intersection before any outcome tests. Record nominal random support, new split-related shortages and shared-reference feasibility in `support_audit.csv`.

- Random arm: sample disjoint cases and references from the full accepted same-category pool without post-hoc pose balancing.
- Imbalanced arm: sample cases from the high half and references from the low half, without replacement within each contrast.
- Prefer one shared category reference and globally disjoint pseudo-cases. If infeasible, use per-pair disjoint draws and explicitly allow/record reuse across pairs.
- Keep original quotas: African 52, Asian 190, European 647. European disjoint placebo is impossible when its original quota exhausts its eligible pool; its 57 cells retain p=1 and are labeled unavailable, never evidence of resistance to pose imbalance.
- Independent deterministic substreams use `SeedSequence([seed, repetition, category_index, arm_index])`. Pairing identifies the shared repetition/design, not identical photo draws between arms.

Every draw must realize positive mean and median absolute-pitch separation in the imbalanced arm. A failure blocks rather than redrawing until a favorable contrast occurs. Pose metrics include pitch quartiles, minimum/range separation, signed pitch, signed/absolute yaw and roll, and available recorded-sex fractions. Category gaps first average cell-specific differences within each repetition and then summarize across repetitions; they are not patient-level estimates.

## Statistics and validation

Use the exact ordered 120 features and original two-sided Mann–Whitney semantics: exact when the smaller group has at most eight observations and that feature has no ties; otherwise asymptotic with tie and continuity correction. Hedges g is pseudo-case minus reference with pooled sample SD and the exact gamma J correction; zero variance leaves g undefined. Each arm/repetition retains one full 9,600-element BH family, with unsupported p=1 and discovery defined as `supported & q<0.05`.

Checkpoints store supported p/g vectors, pose summaries and selection digests, rather than a dense all-repetition full-family tensor. They contain no serialized photographs. The verifier reconstructs all 2,000 BH families with Statsmodels for a final run, regenerates every allocation and pose metric, checks all selection digests, and compares complete supported vectors from selected repetitions against direct SciPy exact/asymptotic tests and an independent `math.lgamma` Hedges calculation. It checks the original BH baseline, result summaries, recurrence, category support, source/code/checkpoint/output hashes, and protected files. Tests use fabricated arrays only and include pose-associated and pose-independent features, ties, zero variance, unsupported cells, deterministic seeds, fallback sampling, corrupted checkpoints and figure rendering.

All tracked files outside this module, relevant prior private baseline/sensitivity/diagnostic/figure outputs, source inputs and existing research/manuscript/publication artifacts are protected by before/after hashes. The module itself is separately hashed for numerical execution and reporting provenance. These entrypoints contain no commit, push or upload action. Only implementation code, this README and fabricated tests may be committed; real-data artifacts remain ignored.

## Private outputs

The output guard requires an untracked, Git-ignored path inside `local_review/`. Every real-data artifact stays there:

- `analysis_plan.json`, `support_audit.csv`, `pitch_availability_audit.csv`: frozen design and feasibility before outcomes.
- `repetition_summary.csv`, `paired_discovery_differences.csv`, `discovery_summary_overall.csv`, `discovery_summary_by_category.csv`, `discovery_summary_by_region.csv`, `discovery_summary_by_cell.csv`: discovery distributions and paired differences. Unique feature-name counts remain distinct from syndrome/feature test counts.
- `pose_gap_summary.csv`, `pose_gap_by_repetition_category.csv`, `repetition_pose_gaps.csv`: aggregate and private cell-specific achieved pose/sex gaps.
- `feature_discovery_frequency.csv`, `test_discovery_frequencies.csv`, `feature_discovery_indicator_correlations.csv`, `pose_discovery_correlations.csv`: descriptive recurrence and correlation audits. Undefined correlations remain missing.
- `pitch_placebo_distribution.{png,pdf,svg}`, `pitch_imbalance_diagnostics.{png,pdf,svg}` and corresponding plot-value CSVs: 300 DPI PNG and vector figures, with European support labeled N/A.
- `RESULTS_AND_LIMITATIONS.md`: full methods/results/limitations, three manuscript-ready paragraphs, exact category distributions, captions and historical benchmark comparison. No actual manuscript file is edited.
- `aggregate_summary.json`, `protected_file_hashes.json`, `validation_report.json`, `source_and_execution_manifest.json`, `checkpoints/`: provenance, checks and reproducible numerical state.

These are sampling-induced same-source pseudo-cohort associations, not clinical syndrome discoveries or causal pitch effects. Other attributes correlated with pose can change simultaneously. Repetitions share a finite image pool and correlated measurements, so simulation quantiles are conditional on this design rather than patient-population confidence intervals. This study gives no empirical p-value for the original GMDB/FairFace contrast and does not establish clinical validity or remove source confounding. The old pooled pitch experiment is discussed only qualitatively; none of its numerical discoveries is reused.
