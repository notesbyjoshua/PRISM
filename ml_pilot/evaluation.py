import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support

from ml_pilot.data_access import DISPLAY_LABELS, LABELS, METADATA, ROOT, PilotError

SCORE_COLUMNS = {"model", "seed", "n", "balanced_accuracy", "accuracy", "macro_f1", "macro_ovr_auroc",
                 "outer_fold_mean_balanced_accuracy", "outer_fold_sd_balanced_accuracy"}
FOLD_COLUMNS = {"model", "seed", "fold", "n_train", "n_test", "inner_seed", "inner_best_balanced_accuracy",
                "selected_parameters", "balanced_accuracy", "accuracy", "macro_f1", "macro_ovr_auroc"}
FOLD_COLUMNS |= {f"{part}_class_{label}_n" for part in ("train", "test") for label in range(5)}
TABLE_COLUMNS = {
    "model_scores.csv": SCORE_COLUMNS, "outer_fold_scores.csv": FOLD_COLUMNS,
    "class_metrics.csv": {"model", "seed", "cohort", "precision", "recall", "f1", "support"},
    "confusion_matrix_counts.csv": {"true_cohort", *DISPLAY_LABELS},
    "confusion_matrix_percent.csv": {"true_cohort", *DISPLAY_LABELS},
    "cohort_summary.csv": {"cohort", "n", "variable", "category", "count", "missing_n", "median", "q25", "q75", "status"},
    "feature_missingness.csv": {"feature", "missing_cells", "n"},
    "coefficient_summary.csv": {"cohort", "feature", "feature_family", "mean_standardized_coefficient",
                                "sd_standardized_coefficient", "sign_consistency_fraction", "n_outer_folds"},
    "family_coefficient_summary.csv": {"cohort", "feature_family", "mean_absolute_standardized_coefficient",
                                       "sd_absolute_standardized_coefficient", "n_features", "n_outer_folds"},
    "repeat_seed_sensitivity.csv": {"model", "n_seeds", "min_fold_mean", "max_fold_mean", "range_fold_means",
                                    "mean_pooled_balanced_accuracy"},
    "permutation_null.csv": {"permutation", "balanced_accuracy"},
}


def restricted_path(path):
    path = Path(path).absolute()
    allowed = ROOT / "local_review"
    if not path.is_relative_to(allowed) or path == allowed:
        raise PilotError("Outputs must be under the ignored repository local_review directory.")
    for ancestor in (path, *path.parents):
        if ancestor == ROOT.parent:
            break
        if ancestor.is_symlink():
            raise PilotError("Symlinks are not allowed in the private output path.")
    if not path.resolve().is_relative_to(allowed.resolve()):
        raise PilotError("Private output path escapes local_review.")
    return path


def prepare_output(path):
    path = restricted_path(path)
    if path.exists():
        raise PilotError("Output directory already exists; choose a fresh run directory to preserve the audit.")
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.mkdir(mode=0o700)
    for name in ("aggregate", "private"):
        (path / name).mkdir(mode=0o700)
    return path


def write_json(value, path):
    path = restricted_path(path)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    path.chmod(0o600)


def write_table(frame, path):
    path = restricted_path(path)
    permitted = TABLE_COLUMNS.get(path.name)
    if not permitted or not set(frame.columns).issubset(permitted) or frame.columns.duplicated().any():
        raise PilotError("Aggregate export rejected: unrecognized table or columns.")
    if any(isinstance(value, (list, dict, tuple, np.ndarray)) for value in frame.to_numpy().ravel()):
        raise PilotError("Aggregate export rejected: nested individual-level arrays are forbidden.")
    frame.to_csv(path, index=False)
    path.chmod(0o600)


def class_tables(y, results):
    records = []
    for result in results:
        precision, recall, f1, support = precision_recall_fscore_support(
            y, result.predictions, labels=range(5), zero_division=0)
        for index, label in enumerate(LABELS):
            records.append({"model": result.model, "seed": result.seed, "cohort": label,
                            "precision": precision[index], "recall": recall[index], "f1": f1[index],
                            "support": int(support[index])})
    return pd.DataFrame(records)


def confusion_tables(y, predictions):
    matrix = confusion_matrix(y, predictions, labels=range(5))
    if matrix.shape != (5, 5) or matrix.sum() != len(y) or (matrix.sum(axis=1) == 0).any():
        raise PilotError("Invalid five-class confusion matrix.")
    counts = pd.DataFrame(matrix, columns=DISPLAY_LABELS)
    percent = pd.DataFrame(100 * matrix / matrix.sum(axis=1, keepdims=True), columns=DISPLAY_LABELS)
    for table in (counts, percent):
        table.insert(0, "true_cohort", DISPLAY_LABELS)
    return counts, percent


def metadata_input(metadata):
    result = metadata.reindex(columns=METADATA).copy()
    for name in METADATA:
        if name == "sex_label":
            result[name] = result[name].where(result[name].isin(["Female", "Male"]), "Missing/other").astype(object)
        else:
            result[name] = pd.to_numeric(result[name], errors="coerce").replace([np.inf, -np.inf], np.nan)
    return result


def cohort_summary(metadata, y, min_cell=10):
    if min_cell < 5:
        raise PilotError("Descriptive cell threshold must be at least five.")
    metadata = metadata_input(metadata)
    records = []
    for label, name in enumerate(LABELS):
        group = metadata.loc[y == label]
        for variable in METADATA:
            base = {"cohort": name, "n": len(group), "variable": variable}
            if variable == "sex_label":
                counts = group[variable].value_counts().reindex(["Female", "Male", "Missing/other"], fill_value=0)
                if ((counts > 0) & (counts < min_cell)).any():
                    records.append({**base, "status": "entire breakdown suppressed: small cell/complement protection"})
                else:
                    records.extend({**base, "category": category, "count": int(count), "status": "aggregate"}
                                   for category, count in counts.items())
            else:
                values = group[variable].dropna()
                missing = len(group) - len(values)
                if len(values) < min_cell or 0 < missing < min_cell:
                    records.append({**base, "status": "summary suppressed: small observed or missing cell"})
                else:
                    records.append({**base, "missing_n": missing, "median": round(float(values.median()), 1),
                                    "q25": round(float(values.quantile(.25)), 1),
                                    "q75": round(float(values.quantile(.75)), 1), "status": "aggregate; rounded"})
    return pd.DataFrame(records)


def coefficient_tables(result, features, definitions):
    coefficients = np.stack(result.coefficients)
    if coefficients.shape[1:] != (5, len(features)):
        raise PilotError("Coefficients do not match the five-class feature schema.")
    families = definitions.set_index("feature").loc[features, "feature_family"]
    records, grouped = [], []
    for label, name in enumerate(LABELS):
        values = coefficients[:, label, :]
        means, stds = values.mean(axis=0), values.std(axis=0, ddof=1)
        consistency = (np.sign(values) == np.sign(means)).mean(axis=0)
        for index, feature in enumerate(features):
            records.append({"cohort": name, "feature": feature, "feature_family": families.iloc[index],
                            "mean_standardized_coefficient": means[index], "sd_standardized_coefficient": stds[index],
                            "sign_consistency_fraction": consistency[index], "n_outer_folds": len(coefficients)})
        for family in families.unique():
            indices = np.flatnonzero(families.to_numpy() == family)
            magnitude = np.abs(values[:, indices]).mean(axis=1)
            grouped.append({"cohort": name, "feature_family": family,
                            "mean_absolute_standardized_coefficient": magnitude.mean(),
                            "sd_absolute_standardized_coefficient": magnitude.std(ddof=1),
                            "n_features": len(indices), "n_outer_folds": len(coefficients)})
    return pd.DataFrame(records), pd.DataFrame(grouped)


def stability_table(scores):
    records = []
    for name, group in scores.groupby("model", sort=False):
        means = group.outer_fold_mean_balanced_accuracy
        records.append({"model": name, "n_seeds": len(group), "min_fold_mean": means.min(),
                        "max_fold_mean": means.max(), "range_fold_means": means.max() - means.min(),
                        "mean_pooled_balanced_accuracy": group.balanced_accuracy.mean()})
    return pd.DataFrame(records)


def figures(y, primary_results, path, permutation):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logistic = next(result for result in primary_results if result.model == "logistic")
    matrix = confusion_matrix(y, logistic.predictions, labels=range(5))
    percentages = matrix / matrix.sum(axis=1, keepdims=True) * 100
    labels = [f"{label}\n(n={n})" for label, n in zip(DISPLAY_LABELS, np.bincount(y))]
    fig, ax = plt.subplots(figsize=(10, 8), layout="constrained")
    plot = ax.imshow(percentages, vmin=0, vmax=100, cmap="Blues")
    for i in range(5):
        for j in range(5):
            ax.text(j, i, f"{matrix[i, j]}\n{percentages[i, j]:.1f}%", ha="center", va="center",
                    color="white" if percentages[i, j] > 55 else "black")
    ax.set(xticks=range(5), yticks=range(5), xticklabels=labels, yticklabels=labels,
           xlabel="Predicted recorded cohort", ylabel="Recorded cohort",
           title=f"Exploratory internal validation: primary logistic model\nOne held-out prediction per patient; n={len(y)}")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.colorbar(plot, ax=ax, label="Within-cohort percentage")
    save_figure(fig, path / "confusion_matrix.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    names = [result.model for result in primary_results]
    scores = [result.score["balanced_accuracy"] for result in primary_results]
    bars = ax.bar(names, scores, color=["#194d75" if name == "logistic" else "#819da3" for name in names])
    for index, result in enumerate(primary_results):
        ax.scatter(np.full(len(result.folds), index), result.folds.balanced_accuracy, color="black", s=18, zorder=3)
    ax.bar_label(bars, fmt="%.3f", padding=4)
    ax.axhline(.2, color="#9e4c35", linestyle="--", label="Constant-class balanced accuracy")
    ax.set(ylim=(0, 1.05), ylabel="Balanced accuracy", xlabel="Declared model; dots are dependent outer-fold scores",
           title=f"Exploratory internal validation; n={len(y)}, five recorded cohorts\n"
                 f"Class sizes: {', '.join(map(str, np.bincount(y)))}; bars use pooled held-out predictions")
    ax.legend(loc="upper left")
    save_figure(fig, path / "model_comparison.png")
    plt.close(fig)
    if permutation["status"] == "EXECUTED":
        fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
        ax.hist(permutation["null_scores"], bins=20, color="#819da3")
        ax.axvline(permutation["observed"], color="#194d75", linewidth=2, label="Observed primary pooled score")
        ax.set(xlabel="Pooled held-out balanced accuracy", ylabel="Label permutations",
               title=f"Exploratory internal permutation null; n={len(y)}\nFull nested tuning repeated under exchangeability")
        ax.legend()
        save_figure(fig, path / "permutation_null.png")
        plt.close(fig)


def save_figure(figure, path):
    path = restricted_path(path)
    figure.savefig(path, dpi=180)
    path.chmod(0o600)


def markdown_table(frame):
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame.itertuples(index=False, name=None):
        cells = [f"{value:.4f}" if isinstance(value, float) else str(value).replace("|", "\\|") for value in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_report(cohort, results, output, manifest, permutation):
    primary_seed = manifest["protocol"]["seed"]
    primary = [result for result in results if result.seed == primary_seed]
    scores = pd.DataFrame([result.score for result in primary])
    classes = class_tables(cohort.y, primary)
    logistic = next(result for result in primary if result.model == "logistic")
    dummy = next(result for result in primary if result.model == "dummy")
    difference = logistic.score["balanced_accuracy"] - dummy.score["balanced_accuracy"]
    interpretation = ("exceeded" if difference > 0 else "equaled" if difference == 0 else "fell below")
    columns = ["model", "balanced_accuracy", "accuracy", "macro_f1", "outer_fold_mean_balanced_accuracy"]
    stability = stability_table(pd.DataFrame([result.score for result in results]))
    permutation_text = "NOT_RUN; no permutation p-value is available."
    if permutation["status"] == "EXECUTED":
        permutation_text = (f"{permutation['n_permutations']} full nested-CV permutations, seed {permutation['seed']}; "
                            f"one-sided Monte Carlo p={permutation['p_value']:.4f}, with the +1 correction. "
                            "This relies on within-dataset exchangeability and does not adjust for unmeasured confounding.")
    text = f"""# PRISM five-cohort ML pilot: exploratory internal validation

## Question and scope
Can 120 fixed FaceKit facial-geometry features distinguish these five recorded GMDB cohort labels? This is a closed-set classification experiment. It cannot screen for disease, identify healthy patients, diagnose unseen syndromes, or establish clinical utility.

## Cohorts
{markdown_table(pd.DataFrame({'exact_recorded_label': LABELS, 'observed_n': np.bincount(cohort.y)}))}

Total n={len(cohort.y)}. The frozen active patient selection was retained: ages 3 <= age < 10, frontal selection, one photograph per nonmissing unique recorded patient ID. There was no new image selection. Source study: {manifest['study_id']}; FaceKit pin: {manifest['facekit_commit']}. Exactly 120 features in the configured schema and S6 order were used. Missing measurement cells: {cohort.audit['missing_feature_cells']}; all preprocessing occurred within training folds. Dataset/schema/code hashes are recorded privately in run_manifest.json.

## Methods
The predeclared primary model was L2 logistic regression (lbfgs). Median imputation, standardization and tuning of C and class weights were contained in the training pipeline. Hyperparameters were selected by inner balanced accuracy, then refitted on the outer training partition. Evaluation used {manifest['protocol']['outer_folds']} shuffled stratified outer folds and {manifest['protocol']['inner_folds']} independently shuffled inner folds, primary outer seed {primary_seed}. Each patient received exactly one held-out prediction per replicate. Patient overlap was prohibited. No global significance rankings, demographic predictors or full-cohort feature selection entered the primary model. Grids and independently derived inner seeds are in the manifest and outer_fold_scores.csv.

The majority-class dummy was evaluated on identical outer splits. SVM/forest, where present, are exploratory comparisons. Metadata-only logistic regression, where present, uses permitted age/recorded sex/pose fields, with imputation, encoding and tuning inside training folds; unknown sex is grouped as Missing/other. No winning model was selected using outer scores.

## Results
{markdown_table(scores[columns])}

The primary pooled held-out balanced accuracy {interpretation} the actually evaluated dummy by {difference:.4f}. Constant-class balanced accuracy is 0.20. Pooled metrics and arithmetic fold means differ slightly because fold class counts differ. Fold dispersion is descriptive, not a confidence interval. Class precision/recall/F1 and support:

{markdown_table(classes.loc[classes.model.eq('logistic'), ['cohort', 'precision', 'recall', 'f1', 'support']])}

Raw and row-normalized confusion matrices are in the aggregate CSVs and confusion_matrix.png. model_comparison.png shows pooled scores and individual fold scores. Macro one-vs-rest AUROC is secondary and appears only where class-aligned predicted probabilities are available. These probabilities were not separately calibrated. SVM decision functions were not treated as probabilities; no SVM AUROC was calculated. Internal AUROC does not establish clinical calibration.

## Stability and confounding
{markdown_table(stability)}

The repeated splits reuse the same patients; their ranges are neither independent replications nor confidence intervals. The default repeat outer seeds are 20261115 and 20261116; the exact executed seeds are in the manifest. Rounded age/sex/pose summaries are in cohort_summary.csv; small cells and potentially revealing complements are suppressed. This conservative suppression is an internal reporting safeguard, not approval for public release. Cohort-level differences and any metadata-only predictability can reflect age or acquisition confounding and do not validate morphology. No source/acquisition-site metadata were added or inferred. Alternate-photo and Noonan-for-MPS sensitivities: NOT_RUN.

Permutation chance check: {permutation_text}

## Interpretation
These estimates concern discrimination among the five included recorded groups only. The confusion matrix describes which included categories were confused; it does not establish the biological cause of those errors. Standardized coefficient means, fold variation and sign consistency are in coefficient_summary.csv, with S6-family magnitude summaries in family_coefficient_summary.csv. Correlated measurements and varying training scales can destabilize coefficients. Family summaries are descriptive model behavior, not independently validated biomarkers or causal anatomical contributions.

## Limitations
The cohort is small (n={len(cohort.y)}), imbalanced and selected for high-count labels. Features are correlated; image-source conditions, age and pose may confound results. Mucopolysaccharidoses is an umbrella recorded cohort category whose subtype definition requires mentor review. There is no independent external cohort, genetic-confirmation audit beyond recorded labels, or diagnostic calibration. Recorded patient identifiers may have unresolved identity-linkage uncertainty. Nested cross-validation reduces tuning leakage but cannot remove dataset-level bias. Comparing exploratory models on these outer results does not yield an unbiased assessment of a retrospectively chosen winner. No confidence interval is claimed from the five dependent fold scores or repeated splits.

## Next research step
Evaluate the prespecified procedure in independent patients, preferably with matched age and acquisition conditions and a broader syndrome differential. Any final all-data fit is a research artifact; its training data are not an independent test set. PI review is required before release of even aggregate results under publication/PI_REVIEW.md. No release was performed.
"""
    path = restricted_path(output / "aggregate/results_report.md")
    path.write_text(text)
    path.chmod(0o600)
