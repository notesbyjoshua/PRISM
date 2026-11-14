from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ml_pilot import cv
from ml_pilot.data_access import PilotError, validate_cohort


def test_patient_splits_cover_once_and_are_reproducible(artificial_cohort):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    splits = cv.patient_splits(cohort.y, cohort.patient_ids, 5, cv.PRIMARY_SEED)
    again = cv.patient_splits(cohort.y, cohort.patient_ids, 5, cv.PRIMARY_SEED)
    coverage = np.zeros(len(cohort.y))
    for (train, test), (a_train, a_test) in zip(splits, again):
        assert not set(cohort.patient_ids[train]) & set(cohort.patient_ids[test])
        np.testing.assert_array_equal(train, a_train)
        np.testing.assert_array_equal(test, a_test)
        coverage[test] += 1
    assert (coverage == 1).all()
    bad = cohort.patient_ids.copy()
    bad[1] = bad[0]
    with pytest.raises(PilotError, match="one row"):
        cv.patient_splits(cohort.y, bad, 5, cv.PRIMARY_SEED)


def test_outer_and_inner_preprocessing_fit_only_training_partitions(artificial_cohort, monkeypatch):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    X = cohort.X.copy()
    X.iloc[0, 0] = np.nan
    splits = cv.patient_splits(cohort.y, cohort.patient_ids, 5, cv.PRIMARY_SEED)
    X.iloc[splits[0][1], 1] += 1e6
    fit_sizes = []
    original = cv.StandardScaler.fit

    def track_fit(self, values, y=None, **kwargs):
        fit_sizes.append(len(values))
        return original(self, values, y, **kwargs)

    monkeypatch.setattr(cv.StandardScaler, "fit", track_fit)
    observed = []

    def inspect_model(estimator, train, test):
        median = np.nanmedian(X.iloc[train], axis=0)
        np.testing.assert_allclose(estimator.named_steps["imputer"].statistics_, median)
        imputed = np.where(np.isnan(X.iloc[train]), median, X.iloc[train])
        np.testing.assert_allclose(estimator.named_steps["scaler"].mean_, imputed.mean(axis=0))
        assert estimator.named_steps["scaler"].n_samples_seen_ == len(train)
        observed.append((train, test))

    result = cv.nested_cv(X, cohort.y, cohort.patient_ids, observer=inspect_model)
    expected_sizes = []
    for fold, (train, test) in enumerate(splits, 1):
        inner = cv.StratifiedKFold(n_splits=3, shuffle=True, random_state=cv.PRIMARY_SEED + 1000 + fold)
        for inner_train, _ in inner.split(X.iloc[train], cohort.y[train]):
            expected_sizes.extend([len(inner_train)] * 8)
        expected_sizes.append(len(train))
    assert sorted(fit_sizes) == sorted(expected_sizes)
    assert max(fit_sizes) < len(X)
    assert len(observed) == 5
    assert len(result.predictions) == len(X)
    assert result.probabilities.shape == (len(X), 5)


def test_seed_reproduces_scores_and_signal_beats_dummy(artificial_cohort):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    a = cv.nested_cv(cohort.X, cohort.y, cohort.patient_ids, outer_folds=2, inner_folds=2)
    b = cv.nested_cv(cohort.X, cohort.y, cohort.patient_ids, outer_folds=2, inner_folds=2)
    dummy = cv.nested_cv(cohort.X, cohort.y, cohort.patient_ids, "dummy", outer_folds=2, inner_folds=2)
    np.testing.assert_array_equal(a.predictions, b.predictions)
    np.testing.assert_allclose(a.probabilities, b.probabilities)
    assert a.score == b.score
    assert a.score["balanced_accuracy"] > dummy.score["balanced_accuracy"] + .5


def test_designed_independent_null_is_near_chance():
    rng = np.random.default_rng(406)
    y = np.repeat(range(5), 100)
    X = pd.DataFrame(rng.normal(size=(len(y), 120)))
    result = cv.nested_cv(X, y, np.arange(len(y)), outer_folds=3, inner_folds=2)
    assert .12 <= result.score["balanced_accuracy"] <= .28


def test_all_missing_training_features_retained():
    estimator, grid = cv.model_spec("logistic", 3)
    y = np.tile(range(5), 5)
    X = np.ones((25, 120))
    X[:, 0] = np.nan
    estimator.fit(X, y)
    assert estimator.named_steps["model"].coef_.shape == (5, 120)


@pytest.mark.parametrize("model", ["svm", "forest", "metadata_logistic"])
def test_comparators_smoke(artificial_cohort, model, monkeypatch):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    original = cv.model_spec

    def small_spec(name, seed):
        estimator, grid = original(name, seed)
        if name == "forest":
            estimator.set_params(model__n_estimators=10)
        return estimator, {key: values[:1] for key, values in grid.items()}

    monkeypatch.setattr(cv, "model_spec", small_spec)
    if model == "metadata_logistic":
        from ml_pilot.evaluation import metadata_input
        X = metadata_input(cohort.metadata)
    else:
        X = cohort.X
    result = cv.nested_cv(X, cohort.y, cohort.patient_ids, model, outer_folds=2, inner_folds=2)
    assert len(result.predictions) == len(cohort.y)
    assert (result.probabilities is None) == (model == "svm")


def test_permutations_rerun_full_procedure_and_preserve_counts(monkeypatch):
    y = np.repeat(range(5), 10)
    seen = []

    def nested_spy(X, labels, ids, name, outer_seed, outer_folds, inner_folds, jobs):
        np.testing.assert_array_equal(np.bincount(labels), np.bincount(y))
        assert name == "logistic" and outer_folds == 5 and inner_folds == 3
        seen.append(labels.copy())
        return SimpleNamespace(score={"balanced_accuracy": .2})

    monkeypatch.setattr(cv, "nested_cv", nested_spy)
    args = (pd.DataFrame(np.zeros((50, 120))), y, np.arange(50), .8)
    result = cv.permutation_null(*args, permutations=100)
    assert len(seen) == 100 and not np.array_equal(seen[0], y)
    assert result["p_value"] == pytest.approx(1 / 101)
    with pytest.raises(PilotError, match="at least 100"):
        cv.permutation_null(*args, permutations=10)
