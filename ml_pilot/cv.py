import json
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import LinearSVC
from threadpoolctl import threadpool_limits

from ml_pilot.data_access import PilotError

PRIMARY_SEED = 20261114
REPEAT_SEEDS = (20261115, 20261116)
C_GRID = [0.01, 0.1, 1, 10]
WEIGHTS = [None, "balanced"]


@dataclass
class CVResult:
    model: str
    seed: int
    predictions: np.ndarray
    probabilities: np.ndarray | None
    folds: pd.DataFrame
    score: dict
    coefficients: list


def preprocessing():
    return [("imputer", SimpleImputer(strategy="median", keep_empty_features=True)), ("scaler", StandardScaler())]


def model_spec(name, seed):
    if name == "dummy":
        return DummyClassifier(strategy="most_frequent", random_state=seed), {}
    if name in ("logistic", "metadata_logistic"):
        model = LogisticRegression(solver="lbfgs", max_iter=5000, tol=1e-6, random_state=seed)
        steps = preprocessing()
        if name == "metadata_logistic":
            numeric = Pipeline(preprocessing())
            categorical = Pipeline([
                ("imputer", SimpleImputer(strategy="constant", fill_value="Missing", keep_empty_features=True)),
                ("encoding", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
            ])
            steps = [("metadata", ColumnTransformer([
                ("numeric", numeric, ["conditional_age_years", "pose_yaw", "pose_pitch", "pose_roll"]),
                ("categorical", categorical, ["sex_label"]),
            ]))]
        return Pipeline(steps + [("model", model)]), {"model__C": C_GRID, "model__class_weight": WEIGHTS}
    if name == "svm":
        model = LinearSVC(dual="auto", max_iter=20000, tol=1e-5, random_state=seed)
        return Pipeline(preprocessing() + [("model", model)]), {"model__C": C_GRID, "model__class_weight": WEIGHTS}
    if name == "forest":
        model = RandomForestClassifier(n_estimators=300, random_state=seed, n_jobs=1)
        grid = {"model__min_samples_leaf": [3, 8], "model__max_features": ["sqrt", 0.5],
                "model__class_weight": WEIGHTS}
        return Pipeline([("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
                         ("model", model)]), grid
    raise PilotError("Unknown declared model.")


def metrics(y, predicted, probabilities=None):
    result = {"balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
              "accuracy": float(accuracy_score(y, predicted)),
              "macro_f1": float(f1_score(y, predicted, labels=np.arange(5), average="macro", zero_division=0))}
    if probabilities is not None:
        values = np.asarray(probabilities)
        if (values.shape != (len(y), 5) or not np.isfinite(values).all() or (values < 0).any()
                or not np.allclose(values.sum(axis=1), 1)):
            raise PilotError("Invalid class-aligned probabilities; AUROC cannot be calculated.")
        result["macro_ovr_auroc"] = float(roc_auc_score(y, values, labels=np.arange(5), multi_class="ovr", average="macro"))
    return result


def patient_splits(y, patient_ids, n_splits, seed):
    y, patient_ids = np.asarray(y), pd.Series(patient_ids, dtype="string")
    if (len(y) != len(patient_ids) or patient_ids.isna().any() or patient_ids.str.strip().eq("").any()
            or patient_ids.str.strip().duplicated().any()):
        raise PilotError("CV requires one row per nonmissing unique recorded patient ID.")
    if set(y) != set(range(5)) or n_splits < 2 or np.bincount(y).min() < n_splits:
        raise PilotError("CV requires exactly five classes and enough patients per class for every fold.")
    splits = list(StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))
    coverage = np.zeros(len(y), dtype=int)
    for train, test in splits:
        if set(patient_ids.iloc[train]) & set(patient_ids.iloc[test]):
            raise PilotError("Recorded patient overlap between training and held-out fold.")
        coverage[test] += 1
    if not np.all(coverage == 1):
        raise PilotError("Outer held-out predictions must cover each recorded patient exactly once.")
    return splits


def tune_model(X, y, name, seed, inner_folds, jobs):
    if inner_folds < 2 or np.bincount(y).min() < inner_folds:
        raise PilotError("Insufficient training patients per class for inner CV.")
    estimator, grid = model_spec(name, seed)
    inner = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=seed)
    with warnings.catch_warnings(), threadpool_limits(limits=1):
        warnings.simplefilter("error", ConvergenceWarning)
        if grid:
            search = GridSearchCV(estimator, grid, scoring="balanced_accuracy", cv=inner, refit=True,
                                  n_jobs=jobs, error_score="raise", return_train_score=False)
            search.fit(X, y)
            return search.best_estimator_, search.best_params_, float(search.best_score_)
        estimator.fit(X, y)
    return estimator, {}, None


def nested_cv(X, y, patient_ids, name="logistic", seed=PRIMARY_SEED, outer_folds=5, inner_folds=3,
              jobs=1, progress=None, observer=None):
    X = X.reset_index(drop=True)
    y = np.asarray(y, dtype=int)
    if len(X) != len(y):
        raise PilotError("Predictor and target row counts differ.")
    splits = patient_splits(y, patient_ids, outer_folds, seed)
    predicted = np.full(len(y), -1, dtype=int)
    coverage = np.zeros(len(y), dtype=int)
    probabilities = None if name in ("svm", "dummy") else np.full((len(y), 5), np.nan)
    rows, coefficients = [], []
    for fold, (train, test) in enumerate(splits, 1):
        inner_seed = seed + 1000 + fold
        estimator, parameters, inner_score = tune_model(X.iloc[train], y[train], name, inner_seed, inner_folds, jobs)
        predicted[test] = estimator.predict(X.iloc[test])
        coverage[test] += 1
        probability = None
        if probabilities is not None:
            probability = estimator.predict_proba(X.iloc[test])[:, np.argsort(estimator.classes_)]
            probabilities[test] = probability
        row = {"model": name, "seed": seed, "fold": fold, "n_train": len(train), "n_test": len(test),
               "inner_seed": inner_seed, "inner_best_balanced_accuracy": inner_score,
               "selected_parameters": json.dumps(parameters, sort_keys=True), **metrics(y[test], predicted[test], probability)}
        for label in range(5):
            row[f"train_class_{label}_n"] = int(np.sum(y[train] == label))
            row[f"test_class_{label}_n"] = int(np.sum(y[test] == label))
        rows.append(row)
        if name == "logistic":
            coefficients.append(estimator.named_steps["model"].coef_[np.argsort(estimator.classes_)].copy())
        if observer:
            observer(estimator, train.copy(), test.copy())
        if progress:
            progress(f"{name}, seed {seed}: completed outer fold {fold}/{outer_folds}")
    if not np.all(coverage == 1) or (predicted < 0).any():
        raise PilotError("Invalid out-of-fold coverage.")
    folds = pd.DataFrame(rows)
    score = {"model": name, "seed": seed, "n": len(y), **metrics(y, predicted, probabilities),
             "outer_fold_mean_balanced_accuracy": float(folds.balanced_accuracy.mean()),
             "outer_fold_sd_balanced_accuracy": float(folds.balanced_accuracy.std(ddof=1))}
    return CVResult(name, seed, predicted, probabilities, folds, score, coefficients)


def permutation_null(X, y, patient_ids, observed, permutations=100, seed=7312026,
                     outer_seed=PRIMARY_SEED, outer_folds=5, inner_folds=3, jobs=1, progress=None):
    if permutations < 100:
        raise PilotError("The declared chance check requires at least 100 full nested-CV label permutations.")
    rng = np.random.default_rng(seed)
    scores = []
    for index in range(permutations):
        shuffled = rng.permutation(y)
        result = nested_cv(X, shuffled, patient_ids, "logistic", outer_seed, outer_folds, inner_folds, jobs)
        scores.append(result.score["balanced_accuracy"])
        if progress:
            progress(f"Completed full nested-CV label permutation {index + 1}/{permutations}")
    p_value = (1 + sum(score >= observed for score in scores)) / (permutations + 1)
    return {"status": "EXECUTED", "n_permutations": permutations, "seed": seed, "observed": observed,
            "p_value": p_value, "null_scores": scores,
            "null_assumption": "within-dataset label exchangeability; full procedure retuned for each permutation"}
