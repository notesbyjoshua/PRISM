# PRISM five-cohort machine-learning pilot

This isolated analysis asks whether 120 existing FaceKit geometry measurements distinguish five recorded GMDB cohort labels. The primary sample is 246 patients: Williams-Beuren (67), Cornelia de Lange (64), Coffin-Siris (43), Mucopolysaccharidoses (37), and Kabuki (35). Exact identifiers and their fixed order are defined in `data_access.py`. Mucopolysaccharidoses is an umbrella recorded category, pending mentor clarification of subtypes.

This is exploratory internal, closed-set evaluation. It does not identify healthy patients, diagnose unfamiliar syndromes, establish clinical calibration, or provide external validation. The existing 48-syndrome significance analysis and manuscript are untouched.

## Environment and checks

Use Python 3.11 in a separate environment, or the existing `.venv` when its versions already match. The publication dependency file is unchanged.

```sh
python3 -m venv /tmp/prism-ml-venv
/tmp/prism-ml-venv/bin/python -m pip install -r ml_pilot/requirements.txt
/tmp/prism-ml-venv/bin/python -m pytest ml_pilot/tests -q
/tmp/prism-ml-venv/bin/python -m ml_pilot.train --dry-run
```

Keep any newly created environment outside version control; the example uses a disposable location outside this checkout. The existing `.venv/` is already ignored. Installation is unnecessary when the existing environment matches the pinned versions.

The dry run checks the active configuration and file availability, prints planned labels/counts, and never reads patient rows, validates observed counts, or fits a model. An absent current private input yields `BLOCKED_PRIVATE_DATA_MISSING`, exit code 2. No aggregate tables, historical 125-feature exports, FairFace rows or fabricated clinical rows are used as substitutes.

## Authorization and private execution

Before patient-level execution, confirm with the PI/mentor that GMDB research authorization and IRB scope cover exploratory syndrome classification. The CLI requires a local JSON record of that confirmation. An author's explicit attestation can be recorded as such; the software cannot independently verify a protocol. Do not invent an approval identifier or approver. Record actual conditions, and keep the record private. This command supports local internal analysis only. Even aggregate release remains subject to `publication/PI_REVIEW.md`.

Create a permission-restricted file under ignored `local_review/`, with this structure and actual approval details:

```json
{
  "analysis": "PRISM_FIVE_SYNDROME_ML_PILOT",
  "approved": true,
  "gmdb_research_authorization_confirmed": true,
  "irb_scope_confirmed": true,
  "approval_reference": "Actual approval reference or explicit author attestation; no invented protocol ID",
  "approved_on": "YYYY-MM-DD",
  "use_conditions": "Actual permitted use conditions",
  "internal_local_analysis_only": true
}
```

Then run:

```sh
.venv/bin/python -m ml_pilot.train --cohort largest5 --seed 20261114 \
  --outer-folds 5 --inner-folds 3 \
  --authorization local_review/ml_pilot_approval.json \
  --output-dir local_review/ml_pilot_20261114
```

This evaluates dummy, primary logistic and exploratory linear SVM models at three outer seeds: 20261114, 20261115, 20261116. Use a new output directory for every run; existing results are never overwritten. Patient independence, exact cohort counts, age eligibility, frontal selection, numeric measurement types and schema/provenance are validated before fitting. Missing feature cells are retained and reported, with median imputation confined to training folds. Infinite values cause a stop. A feature entirely missing in a training partition is retained with zero imputation; it is not selected or dropped using held-out information.

The full declared extension command is:

```sh
.venv/bin/python -m ml_pilot.train --cohort largest5 --seed 20261114 \
  --authorization local_review/ml_pilot_approval.json \
  --forest --metadata-model --permutations 100 --permutation-seed 7312026 --fit-final \
  --output-dir local_review/ml_pilot_20261114_full
```

`--forest` adds a limited random forest comparison. `--metadata-model` adds a separate age/sex/pose model as a confounding diagnostic; the primary predictors remain exactly the 120 measurements. Metadata categories are encoded inside folds. Rounded metadata summaries suppress cells below 10 and suppress a whole sex breakdown if a small count could be derived from complements. This is a conservative internal safeguard, not a substitute for source-specific release rules.

`--permutations` is 0 by default (`NOT_RUN`, no p-value). Values from 1 to 99 are rejected. With 100 or more, labels are shuffled across the dataset, preserving class counts, and the entire logistic tuning/evaluation procedure is rerun for every shuffle. Permuted targets are used both for fitting and evaluation. The outer/inner stratification is reconstructed from the permuted labels using the same seed rules. The test statistic is pooled out-of-fold balanced accuracy; the one-sided Monte Carlo p-value is `(1 + null scores >= observed)/(B + 1)`. This is a within-dataset exchangeability test, not proof against confounding.

`--fit-final` tunes logistic regression with three inner folds on all eligible patients after evaluation and saves a private fitted pipeline. Its inner selection score is not test performance. The final artifact contains the pipeline, feature order and integer-to-exact-label mapping; no patient IDs or prediction arrays are saved. The code never uploads models or registers them externally. Load joblib artifacts only from trusted local runs.

`--no-repeats` omits the two predefined stability seeds. `--jobs 1` is the tested deterministic default; increase only if local compute warrants it. Numerical linear-algebra threads are limited to one during fitting. The real-data primary protocol requires 5 outer and 3 inner folds; smaller folds occur only in synthetic tests. Alternate-photo and Noonan-replacement experiments are not implemented or run and are explicitly marked `NOT_RUN`.

## Prespecified optimization and evaluation

- Primary: multinomial L2 logistic regression using `lbfgs`, tolerance `1e-6`, maximum 5,000 iterations. Scikit-learn 1.9 defaults to `l1_ratio=0.0`; deprecated `penalty` and `multi_class` arguments are omitted.
- Pipeline: training-fold median imputation, standard scaling, classifier. Hyperparameters: `C=[0.01, 0.1, 1, 10]`, `class_weight=[null, balanced]`. Inner selection uses balanced accuracy; the entire pipeline is refitted in each inner split.
- Majority baseline: `DummyClassifier(strategy="most_frequent")`, on identical outer partitions. Its actual held-out scores are calculated. Constant-class balanced accuracy is 0.20.
- Linear SVM: same imputation/scaling and eight-candidate grid, `dual="auto"`, 20,000 iterations. Decision scores are not passed off as probabilities; SVM AUROC is omitted.
- Optional random forest: 300 trees; `min_samples_leaf=[3,8]`, `max_features=[sqrt,0.5]`, `class_weight=[null,balanced]`. Eight declared candidates; this remains exploratory.
- Outer folds: shuffled stratification with specified seed. Inner seed is `outer_seed + 1000 + one_based_outer_fold`. Inner fold count is three. The fixed feature list and transformations are never optimized on outer outcomes.
- Primary metric: balanced accuracy of pooled held-out predictions. Also report accuracy, macro-F1, per-class precision/recall/F1/support, count and row-normalized confusion matrices, fold scores, arithmetic fold means and descriptive fold SD. Fold mean and pooled balanced accuracy can differ because class counts differ between folds.
- Secondary AUROC: macro one-vs-rest from class-aligned probability vectors only, without separate calibration. Logistic/forest/metadata AUROC is secondary and does not establish clinical calibration.
- Three split replicates reuse patients. Report each score and the range of fold means; this range and fold SD are not confidence intervals. The primary confusion matrix uses only the primary seed.

Optimization ends at these predeclared grids. Ranking exploratory models by outer scores and then treating the winner's score as unbiased is not supported. No full-cohort feature selection, SMOTE, CNN, normalization changes, significance-derived features or deep-learning search is performed. This follows the [scikit-learn nested-CV guidance](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html); implementation was checked against the [LogisticRegression API](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html).

## Outputs and provenance

All run outputs are confined to ignored `local_review/`; symlink paths and existing destinations are rejected. New run directories use mode 0700, files 0600. The CLI applies umask 0077. Intermediate predictions/probabilities remain in memory. Default CSV exports use explicit table/column allowlists; no patient IDs, image IDs, raw features, paths or per-person predictions are exported.

The private run manifest records the active study and FaceKit/model provenance, master/schema/S6/config hashes, pilot code hashes, Git HEAD, package versions, authorization attestation, exact grids, seeds, class counts and missingness. Source hashes are checked before/after loading and after fitting; source or code changes invalidate the run. Failed runs are marked `FAILED_DO_NOT_USE_PARTIAL_OUTPUTS`.

`aggregate/` contains model/fold/class scores, confusion CSVs/figure, comparison figure, conservative metadata summaries, feature missingness counts, fold coefficient summaries, S6-family coefficient summaries, repeat-seed stability, and a qualified `results_report.md`. The optional permutation check adds its null CSV/figure and manifest p-value. `private/final_logistic.joblib` exists only with `--fit-final`. All of these remain internal artifacts requiring PI review before any release. Coefficients summarize model behavior; correlated features and training-scale differences preclude interpreting them as independently validated biomarkers.

## Tests

Tests use wholly artificial fixtures. They check exact schema order, exclusion of hidden metadata, label/count/age/frontal validation, duplicate/missing identities, fold coverage, deterministic seeds, training-only outer and inner preprocessing, synthetic signal/null behavior, comparator execution, full-procedure permutation calls, five-class metric ordering, privacy allowlists, permissions, blocked missing inputs and protected study-file integrity. Synthetic performance is only a check of software mechanics and is never reported as clinical evidence.
