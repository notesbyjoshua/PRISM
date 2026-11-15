# Experimental race-stratified PRISM redesign

This study implements the November 15 research task document using current real GMDB and FairFace
120-feature measurements. The direct user request defers sensitivity analysis, so reference replay,
all-reference inferential alternatives, narrower-Asian alternatives, tiny-cell permutation tests,
bootstrap intervals and repeated-CV/permutation analyses are **NOT RUN**. Existing sensitivity code
and outputs are preserved. Counts and covariate availability are descriptive audits.

The replacement main and demographic analyses are reviewable experimental outputs. They do not
promote themselves over `ACTIVE_ANALYSIS.json`, rewrite existing q-values, change validated cohorts,
edit the manuscript or publish data. All outputs, including small-cell aggregates, remain under
verified Git-ignored `local_review/`. No cloud services or image extraction are used.

## Frozen methods

- Use all 1,845 current selected GMDB records, one image per recorded patient, current QC,
  conditional age 3 to <10 and 120 finite measurements. Reconcile the selected-image manifest,
  accepted parquet cohorts and original pinned extraction. Recount all exact syndrome identifiers;
  require >=5 unique patients in a syndrome/category cell, independently of old >=10 total eligibility.
- Use the current 5,864 accepted FairFace age-category `3-9` paths. Before statistical testing,
  exclude all copies of a byte-identical image with conflicting original race annotations, then
  retain the lexically first path for each remaining repeated image hash. This avoids counting
  exact duplicate photographs as new evidence. Distinct content still does not establish distinct people.
- Preserve the old broad crosswalk: GMDB African/Asian/European; FairFace Black→African,
  White→European, East/Southeast/Central/West Asian, Indian and Middle Eastern→Asian, including
  documented spelling variants. Others/unknown and Latino/Hispanic are excluded. Inclusion of
  Middle Eastern/West Asian in Asian is a substantive source-harmonization choice. The narrower
  definition is deferred with the sensitivity work. Labels are recorded annotations, not genetic ancestry.
- Composition denominator: patients in eligible syndrome/category cells only. Reference maximum
  3,409; scan downward for the largest exactly feasible largest-remainder allocation, with alphabetical
  category ties. Sample without replacement from sorted paths with NumPy seed `20261115`. Reuse
  the one global sample's same-category slice in every matching cell. Require reference n>=20.
- Analysis A: all accepted content-screened mapped FairFace images, three-group Kruskal–Wallis,
  one 120-feature BH family; no pairwise post-hoc family. Report means, SDs, medians, IQRs and nulls.
- Primary: every eligible cell × every ordered feature; two-sided Mann–Whitney with per-feature
  SciPy `auto`-equivalent selection (exact when min n<=8 and no ties; otherwise asymptotic with
  tie/continuity corrections). One new global BH family; unsupported references retain p=1.
  Complete finite measurements are mandatory; unexpected missingness aborts source validation.
  Significance q<.05. Hedges' g uses case-minus-reference sign, pooled sample SD and the exact
  gamma-function correction J(df), rather than its usual approximation. Zero pooled variance
  produces a flagged undefined g. The separate |g|>=.8 description is not a clinical threshold.
- ML: rank global top-five cells by descending count, exact syndrome ID, category. Separately
  select five distinct syndromes in one shared category, preferring n>=15 then requiring n>=10;
  choose maximum total n, largest minimum n, category name. Freeze before significance tests/fitting.
  Five outer/three inner stratified folds, one row per patient; 120 features only; all preprocessing
  inside training pipelines. Compare majority dummy, standardized multinomial L2 logistic,
  standardized RBF SVM, and 200-tree random forest using compact grids recorded in the manifest.
  Score balanced accuracy, macro-F1, accuracy, precision/recall, confusion and genuine top-k scores.
  Fold SD/range are descriptive; no independence-based confidence interval. Report every model.

The current conditional GMDB age interpretation is inherited from the validated cohort and has
not been independently established from raw clinical records. FairFace has age bins rather than
comparable exact ages and is not clinically certified healthy. Matching recorded categories cannot
remove the perfect case/source association, referral effects, hidden demographics or capture artifacts.
Ratios can change through their numerators, denominators or acquisition. ML is an exploratory
internal benchmark without independent clinical validation or calibrated diagnostic probabilities.

## Run from repository root

The existing `.venv` already supplies the pinned dependencies below; no existing dependency files
are changed. To use a separate environment, create one under ignored `local_review/` and install
`revision_results/race_stratified_20261115/requirements-ml.txt` there.

```sh
.venv/bin/python revision_results/race_stratified_20261115/run_race_stratified.py --dry-run
.venv/bin/python revision_results/race_stratified_20261115/run_five_syndrome_ml.py --dry-run
.venv/bin/python revision_results/race_stratified_20261115/run_race_stratified.py
.venv/bin/python revision_results/race_stratified_20261115/run_five_syndrome_ml.py
.venv/bin/python -m pytest revision_results/race_stratified_20261115/tests -q \
  --junitxml=local_review/race_stratified_20261115/tests.xml
.venv/bin/python revision_results/race_stratified_20261115/verify_and_report.py
```

Both analysis entrypoints accept `--seed`, `--output-dir` and `--dry-run`; statistics also accepts
`--reference-size` as a maximum. Use a fresh ignored output directory for changed seeds or parameters,
and use the same seed/output for the paired statistics and ML run. Dry-run writes only private source,
count, mapping and selection audits, without patient rows or inferential tests. Missing or invalid
authorized sources produce explicit BLOCKED status; historical and synthetic fallback is forbidden.
The original numerical arrays are validated exactly using round-trip CSV parsing against accepted
parquet and original extraction arrays. No tolerances conceal changed measurements.

Default private outputs include `RESULTS_REPORT.md`, complete primary and Analysis A p/q/g tables,
selection and source audits, formula/unit metadata, aggregate plots, CV summaries, all four 5×5
confusion matrices in long format, descriptive logistic coefficients, and private patient folds and OOF
predictions. `PRIVATE_` files must never be published. The feature matrix contains no metadata and
is also private. Models are not serialized. Public code and fabricated test fixtures are the only
new repository artifacts. The report distinguishes COMPLETED, BLOCKED, NOT APPLICABLE and NOT RUN.

Manifests record input hashes, feature order, code hashes (including uncommitted code), Git revision,
dependencies, seeds, eligibility, exclusions, test denominators, parameters and output hashes.
The verifier reconstructs BH, rechecks the sampling, independently checks stored U statistics,
rebuilds ML metrics/confusion from private OOF predictions, and confirms protected files/source hashes.
The final report and validation summary are generated after execution. PI/privacy review is required
before any aggregate release, especially small cells; these tools contain no publication action.

Numerical method references: [SciPy Mann–Whitney](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
and [scikit-learn nested CV](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html).
