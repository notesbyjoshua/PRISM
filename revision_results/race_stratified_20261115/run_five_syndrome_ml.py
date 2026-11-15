import argparse
import json
import re
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score,
                             precision_recall_fscore_support, top_k_accuracy_score)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from threadpoolctl import threadpool_limits

from common import (PRIVATE, SEED, cohort_counts, load_sources, private_output, protected_hashes, provenance,
                    select_ml, sha, validate_features, validate_patients, write_json)
from run_race_stratified import selection_audit


def patient_splits(y, patients, folds, seed):
    y, patients = np.asarray(y), pd.Series(patients, dtype='string')
    if (len(y) != len(patients) or patients.isna().any() or patients.str.strip().eq('').any()
            or patients.duplicated().any()):
        raise ValueError('CV requires one row per nonmissing recorded patient')
    if set(y) != set(range(5)) or np.bincount(y).min() < folds:
        raise ValueError('Five supported classes required in every fold')
    splits = list(StratifiedKFold(folds, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))
    coverage = np.zeros(len(y), dtype=int)
    for train, test in splits:
        assert not set(patients.iloc[train]) & set(patients.iloc[test])
        assert set(y[train]) == set(y[test]) == set(range(5))
        coverage[test] += 1
    assert np.all(coverage == 1)
    return splits


def model_spec(name, seed):
    if name == 'dummy':
        return DummyClassifier(strategy='most_frequent', random_state=seed), {}
    impute = [('imputer', SimpleImputer(strategy='median', keep_empty_features=True))]
    if name == 'logistic':
        model = LogisticRegression(solver='lbfgs', max_iter=5000, tol=1e-6, random_state=seed)
        steps = impute + [('scaler', StandardScaler()), ('model', model)]
        grid = {'model__C': [.01, .1, 1., 10.], 'model__class_weight': [None, 'balanced']}
    elif name == 'rbf_svm':
        model = SVC(kernel='rbf', decision_function_shape='ovr', break_ties=True, random_state=seed)
        steps = impute + [('scaler', StandardScaler()), ('model', model)]
        grid = {'model__C': [.1, 1., 10.], 'model__gamma': ['scale', .01],
                'model__class_weight': [None, 'balanced']}
    elif name == 'forest':
        model = RandomForestClassifier(n_estimators=200, random_state=seed, n_jobs=1)
        steps = impute + [('model', model)]
        grid = {'model__max_depth': [None, 6], 'model__min_samples_leaf': [1, 3],
                'model__class_weight': [None, 'balanced_subsample']}
    else:
        raise ValueError('Unknown model')
    return Pipeline(steps), grid


def metrics(y, predicted, scores=None):
    result = {'balanced_accuracy': float(balanced_accuracy_score(y, predicted)),
              'macro_f1': float(f1_score(y, predicted, average='macro', labels=range(5), zero_division=0)),
              'accuracy': float(accuracy_score(y, predicted))}
    if scores is not None:
        if scores.shape != (len(y), 5) or not np.isfinite(scores).all():
            raise ValueError('Invalid class-aligned ranking scores')
        for k in [2, 3]:
            result[f'top_{k}_accuracy'] = float(top_k_accuracy_score(y, scores, k=k, labels=np.arange(5)))
    return result


def nested_cv(X, y, patients, name, seed=SEED, observer=None, progress=None):
    y = np.asarray(y, dtype=int)
    outer = patient_splits(y, patients, 5, seed)
    predictions = np.full(len(y), -1, dtype=int)
    scores = None if name == 'dummy' else np.full((len(y), 5), np.nan)
    assignments, fold_rows, inner_rows, coefficients = [], [], [], []
    for fold, (train, test) in enumerate(outer, 1):
        inner_seed = seed + 1000 + fold
        inner = patient_splits(y[train], np.asarray(patients)[train], 3, inner_seed)
        for inner_fold, (inner_train, inner_test) in enumerate(inner, 1):
            assert not set(train[inner_train]) & set(test)
            assert not set(train[inner_test]) & set(test)
            for part, ids in [('train', inner_train), ('validation', inner_test)]:
                inner_rows.extend({'model': name, 'outer_fold': fold, 'inner_fold': inner_fold, 'role': part,
                                   'patient_id': str(np.asarray(patients)[train[i]])} for i in ids)
        estimator, grid = model_spec(name, inner_seed)
        with warnings.catch_warnings(), threadpool_limits(limits=1):
            warnings.simplefilter('error', ConvergenceWarning)
            if grid:
                search = GridSearchCV(estimator, grid, scoring='balanced_accuracy', cv=inner, n_jobs=1,
                                      refit=True, error_score='raise', return_train_score=False)
                search.fit(X.iloc[train], y[train])
                estimator = search.best_estimator_
                parameters, inner_score = search.best_params_, float(search.best_score_)
            else:
                estimator.fit(X.iloc[train], y[train])
                parameters, inner_score = {}, None
            predictions[test] = estimator.predict(X.iloc[test])
            if scores is not None:
                assert np.array_equal(estimator.classes_, np.arange(5))
                values = estimator.decision_function(X.iloc[test]) if name == 'rbf_svm' else estimator.predict_proba(X.iloc[test])
                scores[test] = values
        row = {'model': name, 'fold': fold, 'seed': seed, 'inner_seed': inner_seed,
               'n_train': len(train), 'n_test': len(test), 'inner_best_balanced_accuracy': inner_score,
               'selected_parameters': json.dumps(parameters, sort_keys=True),
               **metrics(y[test], predictions[test], None if scores is None else scores[test])}
        for label in range(5):
            row[f'train_class_{label}_n'] = int(np.sum(y[train] == label))
            row[f'test_class_{label}_n'] = int(np.sum(y[test] == label))
        fold_rows.append(row)
        assignments.extend({'patient_id': str(np.asarray(patients)[i]), 'outer_fold': fold, 'class_id': int(y[i])}
                           for i in test)
        if name == 'logistic':
            coefficients.append(estimator.named_steps['model'].coef_.copy())
        if observer:
            observer(estimator, train, test)
        if progress:
            progress(f'{name}: completed outer fold {fold}/5', flush=True)
    assert (predictions >= 0).all()
    summary = {'model': name, 'n': len(y), 'seed': seed, **metrics(y, predictions, scores)}
    folds = pd.DataFrame(fold_rows)
    summary.update({'fold_balanced_accuracy_mean': float(folds.balanced_accuracy.mean()),
                    'fold_balanced_accuracy_sd_descriptive': float(folds.balanced_accuracy.std(ddof=1)),
                    'fold_balanced_accuracy_min': float(folds.balanced_accuracy.min()),
                    'fold_balanced_accuracy_max': float(folds.balanced_accuracy.max())})
    return summary, predictions, scores, folds, pd.DataFrame(assignments), pd.DataFrame(inner_rows), coefficients


def plot_results(output, summary, matrices, classes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5), layout='constrained')
    ax.bar(summary.model, summary.balanced_accuracy, color=['#a5a5a5', '#28748c', '#de9f49', '#7565a1'])
    ax.axhline(.2, color='#b54143', linestyle='--', label='Five-class balanced chance (20%)')
    ax.set(ylim=(0, 1), ylabel='Out-of-fold balanced accuracy', title='Single-category five-syndrome internal nested CV')
    for i, row in enumerate(summary.itertuples()):
        ax.text(i, row.balanced_accuracy + .02, f'{row.balanced_accuracy:.3f}', ha='center')
    ax.legend(loc='upper left', fontsize=8)
    fig.savefig(output / 'ml_model_comparison.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), layout='constrained')
    short = [re.sub(r'\s+syndrome\b', '', name.split('__')[0].split(';')[0], flags=re.I) for name in classes]
    for ax, (name, matrix) in zip(axes.flat, matrices.items()):
        ax.imshow(matrix, cmap='Blues', vmin=0, vmax=max(m.max() for m in matrices.values()))
        ax.set_xticks(range(5), short, rotation=50, ha='right', fontsize=7)
        ax.set_yticks(range(5), short, fontsize=7)
        ax.set(xlabel='Predicted syndrome', ylabel='Recorded syndrome', title=name)
        for i in range(5):
            for j in range(5):
                ax.text(j, i, str(matrix[i, j]), ha='center', va='center', color='black')
    fig.suptitle('One out-of-fold prediction per patient; all models shown')
    fig.savefig(output / 'ml_confusion_matrix.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Private single-category five-syndrome nested-CV benchmark')
    parser.add_argument('--seed', type=int, default=SEED)
    parser.add_argument('--output-dir', type=str, default=str(PRIVATE))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    output = private_output(args.output_dir)
    before = protected_hashes()
    try:
        disease, _, features, _, audit, sources = load_sources()
    except (FileNotFoundError, ValueError, AssertionError) as error:
        write_json(output / 'ml_status.json', {'status': 'BLOCKED', 'reason': str(error)})
        print(str(error), file=sys.stderr)
        return 2
    counts = cohort_counts(disease)
    selected = selection_audit(output, counts)
    if selected.empty:
        write_json(output / 'ml_status.json', {'status': 'BLOCKED', 'reason': 'No category supports five classes n>=10'})
        return 2
    race, classes = selected.canonical_race.iloc[0], selected.disease.tolist()
    frame = disease.loc[disease.canonical_race.eq(race) & disease.disease.isin(classes)].sort_values('patient_id').copy()
    validate_patients(frame)
    X = validate_features(frame, features).reset_index(drop=True)
    y = frame.disease.map({label: i for i, label in enumerate(classes)}).to_numpy(int)
    patients = frame.patient_id.to_numpy()
    selection = {'race': race, 'classes': classes, 'class_counts': {label: int(frame.disease.eq(label).sum()) for label in classes}}
    manifest = provenance(sources, {'seed': args.seed, 'outer_folds': 5, 'inner_folds': 3, 'n_jobs': 1,
                                    'selection': selection, 'output_dir': str(output), 'dry_run': args.dry_run})
    manifest['features'] = features
    manifest['model_grids'] = {name: model_spec(name, args.seed)[1] for name in ['dummy', 'logistic', 'rbf_svm', 'forest']}
    manifest['source_audit'] = audit
    manifest['protected_before'] = before
    write_json(output / ('ml_dry_run_manifest.json' if args.dry_run else 'ml_manifest.json'), manifest)
    print(json.dumps({'status': 'FEASIBLE', **selection, 'predictors': len(features), 'n': len(frame)}, indent=2), flush=True)
    if args.dry_run:
        assert protected_hashes() == before
        return 0
    frame[['patient_id', 'image_id', 'image_sha256', 'disease', 'canonical_race']].to_csv(
        output / 'PRIVATE_ml_cohort.csv', index=False)
    X.to_csv(output / 'PRIVATE_ml_feature_matrix.csv', index=False)
    summary_rows, per_class, confusion_rows, all_folds, matrices = [], [], [], [], {}
    for name in ['dummy', 'logistic', 'rbf_svm', 'forest']:
        result = nested_cv(X, y, patients, name, args.seed, progress=print)
        summary, predictions, scores, folds, assignments, inner_assignments, coefficients = result
        summary_rows.append(summary)
        all_folds.append(folds)
        assignments.to_csv(output / f'PRIVATE_{name}_outer_folds.csv', index=False)
        inner_assignments.to_csv(output / f'PRIVATE_{name}_inner_folds.csv', index=False)
        oof = pd.DataFrame({'patient_id': patients, 'true_class': y, 'predicted_class': predictions})
        if scores is not None:
            for label in range(5):
                oof[f'ranking_score_{label}'] = scores[:, label]
        oof.to_csv(output / f'PRIVATE_{name}_oof_predictions.csv', index=False)
        precision, recall, f1, support = precision_recall_fscore_support(y, predictions, labels=range(5), zero_division=0)
        for i, label in enumerate(classes):
            per_class.append({'model': name, 'disease': label, 'precision': precision[i], 'recall': recall[i],
                              'f1': f1[i], 'support': int(support[i])})
        matrix = confusion_matrix(y, predictions, labels=range(5))
        matrices[name] = matrix
        for i in range(5):
            for j in range(5):
                confusion_rows.append({'model': name, 'true_disease': classes[i], 'predicted_disease': classes[j],
                                       'count': int(matrix[i, j])})
        if coefficients:
            values = np.asarray(coefficients)
            rows = [{'disease': label, 'feature': feature, 'mean_standardized_coefficient': float(values[:, i, j].mean()),
                     'fold_sd_descriptive': float(values[:, i, j].std(ddof=1))}
                    for i, label in enumerate(classes) for j, feature in enumerate(features)]
            pd.DataFrame(rows).to_csv(output / 'ml_logistic_coefficients_descriptive.csv', index=False)
        print(json.dumps(summary), flush=True)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output / 'ml_cv_summary.csv', index=False)
    pd.DataFrame(per_class).to_csv(output / 'ml_per_class_summary.csv', index=False)
    pd.DataFrame(confusion_rows).to_csv(output / 'ml_confusion_matrix.csv', index=False)
    pd.concat(all_folds, ignore_index=True).to_csv(output / 'ml_fold_metrics.csv', index=False)
    plot_results(output, summary, matrices, classes)
    methods = '''# Five-syndrome ML methods and limitations

This exploratory benchmark uses five distinct exact GMDB syndrome labels within one recorded category.
Classes are frozen by current counts before fitting or statistical feature testing. Predictors are exactly
120 pinned FaceKit measurements. Patient identifiers, image/source/race/demographic fields, pose and labels
are excluded. Each recorded patient and byte-identical source image appears once. Unknown cross-patient
linkage and phenotype/subtype heterogeneity remain limitations.

One stratified five-fold outer CV supplies one held-out prediction per patient. Each outer training set
uses stratified three-fold inner CV to select hyperparameters by balanced accuracy. Every median imputer
and scaler is inside the estimator pipeline and is fitted only on the relevant training fold. There is no
global feature selection. Standardized logistic regression uses the multinomial lbfgs objective and L2
regularization. The RBF SVM uses OVR decision scores for top-k ranking; these are not probabilities.
The forest uses 200 trees and fold-local median imputation. All grids, seeds and dependencies are in
ml_manifest.json; ties follow deterministic GridSearchCV parameter ordering. No estimators are released.

Primary performance is pooled outer out-of-fold balanced accuracy, alongside macro-F1, accuracy, per-class
precision/recall/support, and all four confusion matrices. Top-2/3 use genuine probability or decision
scores for the three fitted models; the dummy has no reported ranking score. The 20% balanced-chance
line and observed majority dummy are both shown. Fold SD/range describe variation; overlapping training
sets do not justify treating folds as independent replicates or using a naive confidence interval.
Comparing several models on these same outer folds does not provide a fresh test set for a selected winner.
Logistic coefficients summarize fold-specific standardized fits, conditional on correlated geometry ratios;
they are descriptive and do not establish causal anatomy. Probability calibration is unassessed.

These are internal CV estimates on 123 selected patients for the default source cohort (the actual n and
classes are recorded in the manifest). Acquisition, demographic, referral and selection effects can remain
within a recorded category. No independent cohort, clinical validation, images, HPO review, permutation
null, repeated-CV sensitivity or diagnostic-utility assessment is included. Sensitivity analysis is deferred
by the user's direct instruction.
'''
    methods = methods.replace('123 selected patients', f'{len(frame)} selected patients')
    (output / 'ml_methods_and_limitations.md').write_text(methods)
    assert protected_hashes() == before, 'Protected analyses changed during ML execution'
    write_json(output / 'ml_status.json', {'status': 'COMPLETED', **selection, 'n': len(frame),
                                          'protected_files_unchanged': True, 'outer_patient_overlap': 0})
    manifest['output_sha256'] = {p.name: sha(p) for p in sorted(output.iterdir())
                                 if p.is_file() and (p.name.startswith('ml_') or p.name.startswith('PRIVATE_'))
                                 and p.name not in ['ml_manifest.json', 'ml_dry_run_manifest.json']}
    write_json(output / 'ml_manifest.json', manifest)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
