"""

To run the program, run this script:

cd /Users/joshua/Documents/PRISM
source .venv/bin/activate
python stats_testing/statistical_analysis.py

"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, kruskal, ttest_ind
from statsmodels.stats.multitest import multipletests


# =========================
# 1. Load data
# =========================

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "stats_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------------
# CONFIGURE YOUR INPUT FILES HERE
# Use an absolute path, or a path relative to the PRISM project directory.
# -------------------------------------------------------------------------
DISEASE_FILE = "data/phenotypes_all_disease.csv"  # put file here
HEALTHY_FILE = ""  # put file here


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare facial phenotype measurements in healthy and disease groups."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        help="Optional combined CSV containing a disease column.",
    )
    parser.add_argument(
        "--healthy-csv",
        type=Path,
        help="CSV containing healthy phenotype measurements.",
    )
    parser.add_argument(
        "--disease-csv",
        type=Path,
        help="CSV containing disease phenotype measurements.",
    )
    parser.add_argument("--bootstrap-iterations", type=int, default=500)
    parser.add_argument("--rank-bootstrap-iterations", type=int, default=200)
    parser.add_argument("--bootstrap-seed", type=int, default=2026)
    parser.add_argument("--max-missing-fraction", type=float, default=0.30)
    parser.add_argument("--min-group-n", type=int, default=10)
    parser.add_argument("--effect-threshold", type=float, default=0.50)
    parser.add_argument("--kw-effect-threshold", type=float, default=0.06)
    parser.add_argument("--mad-threshold", type=float, default=5.0)
    parser.add_argument("--top-k", type=int, default=10)
    return parser.parse_args()


def resolve_file_path(file_path):
    """Resolve configured relative paths from the PRISM project directory."""
    path = Path(file_path).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path


def load_data(args):
    # Command-line paths override the easy-to-edit constants above. If no
    # command-line paths are supplied, running from an IDE uses those constants.
    healthy_path = args.healthy_csv
    disease_path = args.disease_csv

    if healthy_path is None and disease_path is None:
        if bool(HEALTHY_FILE) != bool(DISEASE_FILE):
            raise ValueError(
                "Set both HEALTHY_FILE and DISEASE_FILE in the Load data section."
            )
        if HEALTHY_FILE and DISEASE_FILE:
            healthy_path = resolve_file_path(HEALTHY_FILE)
            disease_path = resolve_file_path(DISEASE_FILE)

    using_split_files = healthy_path is not None or disease_path is not None

    if using_split_files:
        if healthy_path is None or disease_path is None:
            raise ValueError(
                "Provide both --healthy-csv and --disease-csv, or configure both "
                "file constants in the Load data section."
            )

        healthy_path = resolve_file_path(healthy_path)
        disease_path = resolve_file_path(disease_path)
        for label, path in [("healthy", healthy_path), ("disease", disease_path)]:
            if not path.is_file():
                raise FileNotFoundError(f"The configured {label} CSV does not exist: {path}")

        print(f"Loading healthy data from: {healthy_path}")
        print(f"Loading disease data from: {disease_path}")
        healthy_data = pd.read_csv(healthy_path)
        disease_data = pd.read_csv(disease_path)

        # The source file determines which observations are healthy, so overwrite
        # any existing label in the healthy file.
        healthy_data["disease"] = "healthy"

        # Preserve specific disease names when supplied. Otherwise, analyze all
        # rows from the disease file as one disease group.
        if "disease" not in disease_data.columns:
            disease_data["disease"] = "disease"

        return pd.concat([healthy_data, disease_data], ignore_index=True, sort=False)

    if args.input_csv is not None:
        input_path = resolve_file_path(args.input_csv)
        if not input_path.is_file():
            raise FileNotFoundError(f"The combined input CSV does not exist: {input_path}")
        data = pd.read_csv(input_path)
        if "disease" not in data.columns:
            raise ValueError("The combined input CSV must contain a 'disease' column.")
        return data

    raise ValueError(
        "No input files configured. Set DISEASE_FILE and HEALTHY_FILE near the "
        "top of the Load data section, then run the script again."
    )


args = parse_args()
df = load_data(args)

if args.bootstrap_iterations < 0 or args.rank_bootstrap_iterations < 0:
    raise ValueError("Bootstrap iteration counts cannot be negative.")
if not 0 <= args.max_missing_fraction <= 1:
    raise ValueError("--max-missing-fraction must be between 0 and 1.")
if args.min_group_n < 3:
    raise ValueError("--min-group-n must be at least 3.")

# Keep only frontal images
df = df[df["frontal_ok"] == True].copy()

# Optional: exclude pose variables because they are acquisition/quality features,
# not true facial phenotype features.
DROP_POSE = True

exclude_cols = ["disease", "image_id", "frontal_ok"]

if DROP_POSE:
    exclude_cols += ["pose_yaw", "pose_pitch", "pose_roll"]

feature_cols = [
    c for c in df.columns
    if c not in exclude_cols and pd.api.types.is_numeric_dtype(df[c])
]

print("Number of rows:", len(df))
print("Number of features:", len(feature_cols))
print(df["disease"].value_counts())


# =========================
# 2. Define helper functions
# =========================

def cohens_d(x, y):
    """
    Cohen's d for two independent groups.
    Positive value means x has a larger mean than y.
    """
    x = np.asarray(x.dropna(), dtype=float)
    y = np.asarray(y.dropna(), dtype=float)

    nx = len(x)
    ny = len(y)

    if nx < 2 or ny < 2:
        return np.nan

    sx = np.var(x, ddof=1)
    sy = np.var(y, ddof=1)

    pooled_sd = np.sqrt(((nx - 1) * sx + (ny - 1) * sy) / (nx + ny - 2))

    if pooled_sd == 0:
        return np.nan

    return (np.mean(x) - np.mean(y)) / pooled_sd


def hedges_g(x, y):
    """
    Small-sample corrected Cohen's d.
    Usually better than plain Cohen's d for small disease groups.
    """
    d = cohens_d(x, y)

    x = np.asarray(x.dropna(), dtype=float)
    y = np.asarray(y.dropna(), dtype=float)

    nx = len(x)
    ny = len(y)

    if np.isnan(d) or nx + ny <= 3:
        return np.nan

    correction = 1 - (3 / (4 * (nx + ny) - 9))
    return d * correction


def cliffs_delta(x, y):
    """
    Nonparametric effect size.
    Range: -1 to +1.
    Positive means x tends to have larger values than y.
    """
    x = np.asarray(x.dropna(), dtype=float)
    y = np.asarray(y.dropna(), dtype=float)

    if len(x) == 0 or len(y) == 0:
        return np.nan

    greater = 0
    lesser = 0

    for xi in x:
        greater += np.sum(xi > y)
        lesser += np.sum(xi < y)

    return (greater - lesser) / (len(x) * len(y))


def rank_biserial_from_u(u_stat, n1, n2):
    """
    Rank-biserial correlation from Mann-Whitney U.
    Range: -1 to +1.
    """
    if n1 == 0 or n2 == 0:
        return np.nan

    return (2 * u_stat / (n1 * n2)) - 1


def auc_from_u(u_stat, n1, n2):
    """Directional and direction-agnostic ROC AUC from Mann-Whitney U."""
    if n1 == 0 or n2 == 0 or not np.isfinite(u_stat):
        return np.nan, np.nan
    directional_auc = u_stat / (n1 * n2)
    discrimination_auc = max(directional_auc, 1 - directional_auc)
    return directional_auc, discrimination_auc


def epsilon_squared_kruskal(h_stat, n, k):
    """
    Effect size for Kruskal-Wallis.
    Values closer to 1 mean stronger group differences.
    """
    if n <= k:
        return np.nan

    return (h_stat - k + 1) / (n - k)


def remove_mad_outliers(values, threshold):
    """Remove values farther than ``threshold`` scaled MADs from the median."""
    values = pd.Series(values).dropna().astype(float)
    if values.empty:
        return values
    median = values.median()
    mad = np.median(np.abs(values - median))
    if mad == 0 or not np.isfinite(mad):
        return values
    robust_z = 0.67448975 * np.abs(values - median) / mad
    return values[robust_z <= threshold]


def cliffs_delta_fast(x, y):
    """Cliff's delta without constructing the full pairwise comparison matrix."""
    x = np.asarray(x, dtype=float)
    y = np.sort(np.asarray(y, dtype=float))
    if len(x) == 0 or len(y) == 0:
        return np.nan
    lesser = np.searchsorted(y, x, side="left").sum()
    greater = (len(y) - np.searchsorted(y, x, side="right")).sum()
    return (lesser - greater) / (len(x) * len(y))


def bootstrap_effect_ci(x, y, iterations, rng):
    """Percentile 95% CIs for Hedges' g and Cliff's delta."""
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if iterations == 0 or len(x) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan, np.nan

    g_samples = np.empty(iterations)
    delta_samples = np.empty(iterations)
    for i in range(iterations):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        g_samples[i] = hedges_g(pd.Series(xb), pd.Series(yb))
        delta_samples[i] = cliffs_delta_fast(xb, yb)

    g_valid = g_samples[np.isfinite(g_samples)]
    delta_valid = delta_samples[np.isfinite(delta_samples)]
    g_ci = np.quantile(g_valid, [0.025, 0.975]) if len(g_valid) else [np.nan, np.nan]
    d_ci = (
        np.quantile(delta_valid, [0.025, 0.975])
        if len(delta_valid)
        else [np.nan, np.nan]
    )
    return g_ci[0], g_ci[1], d_ci[0], d_ci[1]


def add_fdr(p_values):
    """Benjamini-Hochberg correction that safely preserves missing p-values."""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p_values), np.nan)
    valid = np.isfinite(p_values)
    if valid.any():
        adjusted[valid] = multipletests(p_values[valid], method="fdr_bh")[1]
    return adjusted


def rank_stability(group_df, reference_df, features, iterations, top_k, rng):
    """Frequency with which features rank in the top K by absolute Hedges' g."""
    if iterations == 0 or not features:
        return {feature: np.nan for feature in features}

    x = group_df[features].to_numpy(dtype=float)
    y = reference_df[features].to_numpy(dtype=float)
    counts = np.zeros(len(features), dtype=int)
    k = min(top_k, len(features))

    for _ in range(iterations):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        nx = np.sum(np.isfinite(xb), axis=0)
        ny = np.sum(np.isfinite(yb), axis=0)
        mean_x = np.nanmean(xb, axis=0)
        mean_y = np.nanmean(yb, axis=0)
        var_x = np.nanvar(xb, axis=0, ddof=1)
        var_y = np.nanvar(yb, axis=0, ddof=1)
        denominator = nx + ny - 2
        pooled = np.sqrt(((nx - 1) * var_x + (ny - 1) * var_y) / denominator)
        correction = 1 - (3 / (4 * (nx + ny) - 9))
        scores = np.abs(((mean_x - mean_y) / pooled) * correction)
        scores[(nx < 2) | (ny < 2) | ~np.isfinite(scores)] = -np.inf
        top_indices = np.argsort(scores)[-k:]
        counts[top_indices[scores[top_indices] > -np.inf]] += 1

    return dict(zip(features, counts / iterations))


def robust_two_group_stats(disease_values, healthy_values, threshold):
    disease_clean = remove_mad_outliers(disease_values, threshold)
    healthy_clean = remove_mad_outliers(healthy_values, threshold)
    if len(disease_clean) < 3 or len(healthy_clean) < 3:
        return {
            "robust_n_disease": len(disease_clean),
            "robust_n_healthy": len(healthy_clean),
            "robust_mannwhitney_u": np.nan,
            "robust_mannwhitney_p": np.nan,
            "robust_roc_auc_directional": np.nan,
            "robust_roc_auc": np.nan,
            "robust_cliffs_delta": np.nan,
            "robust_hedges_g": np.nan,
        }
    u_stat, p_value = mannwhitneyu(
        disease_clean, healthy_clean, alternative="two-sided"
    )
    directional_auc, discrimination_auc = auc_from_u(
        u_stat, len(disease_clean), len(healthy_clean)
    )
    return {
        "robust_n_disease": len(disease_clean),
        "robust_n_healthy": len(healthy_clean),
        "robust_mannwhitney_u": u_stat,
        "robust_mannwhitney_p": p_value,
        "robust_roc_auc_directional": directional_auc,
        "robust_roc_auc": discrimination_auc,
        "robust_cliffs_delta": cliffs_delta_fast(disease_clean, healthy_clean),
        "robust_hedges_g": hedges_g(disease_clean, healthy_clean),
    }


# =========================
# 3. Test 1:
#    Healthy vs all disease
# =========================

binary_results = []

healthy_df = df[df["disease"] == "healthy"]
disease_df = df[df["disease"] != "healthy"]
rng = np.random.default_rng(args.bootstrap_seed)

if healthy_df.empty or disease_df.empty:
    raise ValueError("The input must contain both healthy and disease rows.")

# Save reliability/missingness for every feature in every source group. These
# rows are retained even when a feature is excluded from hypothesis testing.
reliability_rows = []
for label, group_df in [("healthy", healthy_df), ("all_disease", disease_df)]:
    for feature in feature_cols:
        n_valid = int(group_df[feature].notna().sum())
        reliability_rows.append({
            "group": label,
            "feature": feature,
            "n_rows": len(group_df),
            "n_valid": n_valid,
            "n_missing": len(group_df) - n_valid,
            "fraction_missing": 1 - (n_valid / len(group_df)),
            "passes_missingness": 1 - (n_valid / len(group_df)) <= args.max_missing_fraction,
            "passes_min_n": n_valid >= args.min_group_n,
        })

disease_labels = sorted([x for x in df["disease"].dropna().unique() if x != "healthy"])
for disease in disease_labels:
    group_df = df[df["disease"] == disease]
    for feature in feature_cols:
        n_valid = int(group_df[feature].notna().sum())
        reliability_rows.append({
            "group": disease,
            "feature": feature,
            "n_rows": len(group_df),
            "n_valid": n_valid,
            "n_missing": len(group_df) - n_valid,
            "fraction_missing": 1 - (n_valid / len(group_df)),
            "passes_missingness": 1 - (n_valid / len(group_df)) <= args.max_missing_fraction,
            "passes_min_n": n_valid >= args.min_group_n,
        })

reliability_results = pd.DataFrame(reliability_rows)
reliability_results["eligible"] = (
    reliability_results["passes_missingness"] & reliability_results["passes_min_n"]
)
reliability_results.to_csv(OUTPUT_DIR / "feature_reliability.csv", index=False)

for feature in feature_cols:
    healthy_values = healthy_df[feature].dropna()
    disease_values = disease_df[feature].dropna()

    missing_healthy = healthy_df[feature].isna().mean()
    missing_disease = disease_df[feature].isna().mean()
    if (
        len(healthy_values) < args.min_group_n
        or len(disease_values) < args.min_group_n
        or missing_healthy > args.max_missing_fraction
        or missing_disease > args.max_missing_fraction
    ):
        continue

    # Mann-Whitney U test
    try:
        u_stat, mw_p = mannwhitneyu(
            disease_values,
            healthy_values,
            alternative="two-sided",
        )
    except ValueError:
        u_stat, mw_p = np.nan, np.nan

    # Welch t-test, included as a secondary parametric check
    try:
        t_stat, t_p = ttest_ind(
            disease_values,
            healthy_values,
            equal_var=False,
            nan_policy="omit",
        )
    except ValueError:
        t_stat, t_p = np.nan, np.nan

    n_disease = len(disease_values)
    n_healthy = len(healthy_values)
    directional_auc, discrimination_auc = auc_from_u(u_stat, n_disease, n_healthy)

    g_low, g_high, delta_low, delta_high = bootstrap_effect_ci(
        disease_values, healthy_values, args.bootstrap_iterations, rng
    )
    result = {
        "feature": feature,
        "comparison": "all_disease_vs_healthy",
        "n_disease": n_disease,
        "n_healthy": n_healthy,
        "fraction_missing_disease": missing_disease,
        "fraction_missing_healthy": missing_healthy,
        "mean_disease": disease_values.mean(),
        "mean_healthy": healthy_values.mean(),
        "median_disease": disease_values.median(),
        "median_healthy": healthy_values.median(),
        "mean_difference_disease_minus_healthy": disease_values.mean() - healthy_values.mean(),
        "mannwhitney_u": u_stat,
        "mannwhitney_p": mw_p,
        "roc_auc_directional": directional_auc,
        "roc_auc": discrimination_auc,
        "rank_biserial": rank_biserial_from_u(u_stat, n_disease, n_healthy),
        "cliffs_delta": cliffs_delta(disease_values, healthy_values),
        "cliffs_delta_ci_95_low": delta_low,
        "cliffs_delta_ci_95_high": delta_high,
        "cohens_d": cohens_d(disease_values, healthy_values),
        "hedges_g": hedges_g(disease_values, healthy_values),
        "hedges_g_ci_95_low": g_low,
        "hedges_g_ci_95_high": g_high,
        "welch_t": t_stat,
        "welch_p": t_p,
    }
    result.update(robust_two_group_stats(disease_values, healthy_values, args.mad_threshold))
    binary_results.append(result)

binary_results = pd.DataFrame(binary_results)

# FDR correction
binary_results["mannwhitney_q_fdr"] = add_fdr(binary_results["mannwhitney_p"])

binary_results["welch_q_fdr"] = add_fdr(binary_results["welch_p"])
binary_results["robust_mannwhitney_q_fdr"] = add_fdr(
    binary_results["robust_mannwhitney_p"]
)
binary_results["high_priority"] = (
    (binary_results["mannwhitney_q_fdr"] < 0.05)
    & (binary_results["hedges_g"].abs() >= args.effect_threshold)
)

binary_features = binary_results["feature"].tolist()
binary_stability = rank_stability(
    disease_df,
    healthy_df,
    binary_features,
    args.rank_bootstrap_iterations,
    args.top_k,
    rng,
)
binary_results["top_k_bootstrap_frequency"] = binary_results["feature"].map(binary_stability)
binary_results["top_k"] = args.top_k
binary_results["rank_bootstrap_iterations"] = args.rank_bootstrap_iterations

binary_results = binary_results.sort_values(
    ["mannwhitney_q_fdr", "hedges_g"],
    ascending=[True, False],
)

binary_results.to_csv(OUTPUT_DIR / "stats_healthy_vs_disease.csv", index=False)

print("\nTop healthy vs disease features:")
print(binary_results.head(20))


# =========================
# 4. Test 2:
#    Each disease vs healthy
# =========================

pairwise_results = []

for disease in disease_labels:
    disease_specific_df = df[df["disease"] == disease]

    for feature in feature_cols:
        healthy_values = healthy_df[feature].dropna()
        disease_values = disease_specific_df[feature].dropna()

        missing_healthy = healthy_df[feature].isna().mean()
        missing_disease = disease_specific_df[feature].isna().mean()
        if (
            len(healthy_values) < args.min_group_n
            or len(disease_values) < args.min_group_n
            or missing_healthy > args.max_missing_fraction
            or missing_disease > args.max_missing_fraction
        ):
            continue

        try:
            u_stat, mw_p = mannwhitneyu(
                disease_values,
                healthy_values,
                alternative="two-sided",
            )
        except ValueError:
            u_stat, mw_p = np.nan, np.nan

        directional_auc, discrimination_auc = auc_from_u(
            u_stat, len(disease_values), len(healthy_values)
        )

        g_low, g_high, delta_low, delta_high = bootstrap_effect_ci(
            disease_values, healthy_values, args.bootstrap_iterations, rng
        )
        result = {
            "feature": feature,
            "disease": disease,
            "comparison": f"{disease}_vs_healthy",
            "n_disease": len(disease_values),
            "n_healthy": len(healthy_values),
            "fraction_missing_disease": missing_disease,
            "fraction_missing_healthy": missing_healthy,
            "mean_disease": disease_values.mean(),
            "mean_healthy": healthy_values.mean(),
            "median_disease": disease_values.median(),
            "median_healthy": healthy_values.median(),
            "mean_difference_disease_minus_healthy": disease_values.mean() - healthy_values.mean(),
            "mannwhitney_u": u_stat,
            "mannwhitney_p": mw_p,
            "roc_auc_directional": directional_auc,
            "roc_auc": discrimination_auc,
            "rank_biserial": rank_biserial_from_u(
                u_stat,
                len(disease_values),
                len(healthy_values),
            ),
            "cliffs_delta": cliffs_delta(disease_values, healthy_values),
            "cliffs_delta_ci_95_low": delta_low,
            "cliffs_delta_ci_95_high": delta_high,
            "cohens_d": cohens_d(disease_values, healthy_values),
            "hedges_g": hedges_g(disease_values, healthy_values),
            "hedges_g_ci_95_low": g_low,
            "hedges_g_ci_95_high": g_high,
        }
        result.update(
            robust_two_group_stats(disease_values, healthy_values, args.mad_threshold)
        )
        pairwise_results.append(result)

pairwise_results = pd.DataFrame(pairwise_results)

# FDR correction across all disease-vs-healthy feature tests
pairwise_results["mannwhitney_q_fdr_global"] = add_fdr(
    pairwise_results["mannwhitney_p"]
)
pairwise_results["robust_mannwhitney_q_fdr_global"] = add_fdr(
    pairwise_results["robust_mannwhitney_p"]
)

# Also do FDR correction within each disease
pairwise_results["mannwhitney_q_fdr_within_disease"] = np.nan

for disease in disease_labels:
    mask = pairwise_results["disease"] == disease
    if mask.any():
        pairwise_results.loc[mask, "mannwhitney_q_fdr_within_disease"] = add_fdr(
            pairwise_results.loc[mask, "mannwhitney_p"]
        )

pairwise_results["high_priority"] = (
    (pairwise_results["mannwhitney_q_fdr_global"] < 0.05)
    & (pairwise_results["hedges_g"].abs() >= args.effect_threshold)
)
pairwise_results["top_k_bootstrap_frequency"] = np.nan
for disease in disease_labels:
    mask = pairwise_results["disease"] == disease
    if not mask.any():
        continue
    features = pairwise_results.loc[mask, "feature"].tolist()
    stability = rank_stability(
        df[df["disease"] == disease],
        healthy_df,
        features,
        args.rank_bootstrap_iterations,
        args.top_k,
        rng,
    )
    pairwise_results.loc[mask, "top_k_bootstrap_frequency"] = (
        pairwise_results.loc[mask, "feature"].map(stability)
    )
pairwise_results["top_k"] = args.top_k
pairwise_results["rank_bootstrap_iterations"] = args.rank_bootstrap_iterations

# A feature is more disease-specific when it has a qualifying effect in a small
# fraction of the diseases for which it was reliably tested. Specificity alone
# is 1.0 for a feature with no effects, so the weighted score also incorporates
# the median qualifying effect magnitude and becomes zero when none qualify.
specificity_rows = []
for feature, feature_results in pairwise_results.groupby("feature", sort=True):
    qualifying = feature_results[feature_results["high_priority"]]
    diseases_tested = int(feature_results["disease"].nunique())
    diseases_with_effect = int(qualifying["disease"].nunique())
    specificity = (
        1 - diseases_with_effect / diseases_tested if diseases_tested else np.nan
    )
    median_effect = (
        qualifying["hedges_g"].abs().median() if not qualifying.empty else 0.0
    )
    specificity_rows.append({
        "feature": feature,
        "n_diseases_tested": diseases_tested,
        "n_diseases_with_significant_moderate_or_large_effect": diseases_with_effect,
        "n_diseases_with_positive_effect": int((qualifying["hedges_g"] > 0).sum()),
        "n_diseases_with_negative_effect": int((qualifying["hedges_g"] < 0).sum()),
        "fraction_diseases_with_effect": (
            diseases_with_effect / diseases_tested if diseases_tested else np.nan
        ),
        "feature_specificity_score": specificity,
        "median_abs_hedges_g_among_significant": median_effect,
        "max_abs_hedges_g_among_significant": (
            qualifying["hedges_g"].abs().max() if not qualifying.empty else 0.0
        ),
        "specificity_weighted_effect_score": specificity * median_effect,
    })

feature_specificity_results = pd.DataFrame(specificity_rows).sort_values(
    [
        "specificity_weighted_effect_score",
        "n_diseases_with_significant_moderate_or_large_effect",
    ],
    ascending=[False, True],
)
pairwise_results = pairwise_results.merge(
    feature_specificity_results[
        [
            "feature",
            "n_diseases_tested",
            "n_diseases_with_significant_moderate_or_large_effect",
            "fraction_diseases_with_effect",
            "feature_specificity_score",
            "specificity_weighted_effect_score",
        ]
    ],
    on="feature",
    how="left",
)

pairwise_results = pairwise_results.sort_values(
    ["disease", "mannwhitney_q_fdr_within_disease", "hedges_g"],
    ascending=[True, True, False],
)

pairwise_results.to_csv(OUTPUT_DIR / "stats_each_disease_vs_healthy.csv", index=False)

print("\nTop disease-specific features:")
print(pairwise_results.head(30))


# =========================
# 5. Test 3:
#    Kruskal-Wallis across all 6 classes
# =========================

kw_results = []

all_labels = sorted(
    label
    for label, count in df["disease"].value_counts().items()
    if count >= args.min_group_n
)

for feature in feature_cols:
    groups = []
    feature_is_reliable = True

    for label in all_labels:
        label_df = df[df["disease"] == label]
        values = label_df[feature].dropna()
        if (
            len(values) < args.min_group_n
            or label_df[feature].isna().mean() > args.max_missing_fraction
        ):
            feature_is_reliable = False
            break
        groups.append(values)

    if not feature_is_reliable or len(groups) < 2:
        continue

    try:
        h_stat, kw_p = kruskal(*groups)
    except ValueError:
        h_stat, kw_p = np.nan, np.nan

    n_total = sum(len(g) for g in groups)
    k_groups = len(groups)

    robust_groups = [remove_mad_outliers(group, args.mad_threshold) for group in groups]
    if all(len(group) >= 3 for group in robust_groups):
        try:
            robust_h, robust_p = kruskal(*robust_groups)
        except ValueError:
            robust_h, robust_p = np.nan, np.nan
    else:
        robust_h, robust_p = np.nan, np.nan
    robust_n = sum(len(group) for group in robust_groups)

    kw_results.append({
        "feature": feature,
        "n_total": n_total,
        "k_groups": k_groups,
        "kruskal_h": h_stat,
        "kruskal_p": kw_p,
        "epsilon_squared": epsilon_squared_kruskal(h_stat, n_total, k_groups),
        "robust_n_total": robust_n,
        "robust_kruskal_h": robust_h,
        "robust_kruskal_p": robust_p,
        "robust_epsilon_squared": epsilon_squared_kruskal(
            robust_h, robust_n, k_groups
        ),
    })

kw_results = pd.DataFrame(kw_results)

kw_results["kruskal_q_fdr"] = add_fdr(kw_results["kruskal_p"])
kw_results["robust_kruskal_q_fdr"] = add_fdr(kw_results["robust_kruskal_p"])
kw_results["high_priority"] = (
    (kw_results["kruskal_q_fdr"] < 0.05)
    & (kw_results["epsilon_squared"] >= args.kw_effect_threshold)
)

kw_results = kw_results.sort_values(
    ["kruskal_q_fdr", "epsilon_squared"],
    ascending=[True, False],
)

kw_results.to_csv(OUTPUT_DIR / "stats_kruskal_all_groups.csv", index=False)

print("\nTop Kruskal-Wallis features across all groups:")
print(kw_results.head(20))


# =========================
# 6. Optional:
#    Save top significant features
# =========================

SIGNIFICANCE_Q = 0.05

top_binary = binary_results[
    binary_results["mannwhitney_q_fdr"] < SIGNIFICANCE_Q
].copy()

top_pairwise = pairwise_results[
    pairwise_results["mannwhitney_q_fdr_global"] < SIGNIFICANCE_Q
].copy()

top_kw = kw_results[
    kw_results["kruskal_q_fdr"] < SIGNIFICANCE_Q
].copy()

high_priority_binary = binary_results[binary_results["high_priority"]].copy()
high_priority_pairwise = pairwise_results[pairwise_results["high_priority"]].copy()
high_priority_kw = kw_results[kw_results["high_priority"]].copy()

top_binary.to_csv(OUTPUT_DIR / "significant_healthy_vs_disease.csv", index=False)
top_pairwise.to_csv(OUTPUT_DIR / "significant_each_disease_vs_healthy.csv", index=False)
top_kw.to_csv(OUTPUT_DIR / "significant_kruskal_all_groups.csv", index=False)
high_priority_binary.to_csv(
    OUTPUT_DIR / "high_priority_healthy_vs_disease.csv", index=False
)
high_priority_pairwise.to_csv(
    OUTPUT_DIR / "high_priority_each_disease_vs_healthy.csv", index=False
)
high_priority_kw.to_csv(
    OUTPUT_DIR / "high_priority_kruskal_all_groups.csv", index=False
)

rank_stability_results = pd.concat(
    [
        binary_results[["comparison", "feature", "top_k_bootstrap_frequency"]],
        pairwise_results[["comparison", "feature", "top_k_bootstrap_frequency"]],
    ],
    ignore_index=True,
)
rank_stability_results["top_k"] = args.top_k
rank_stability_results["bootstrap_iterations"] = args.rank_bootstrap_iterations
rank_stability_results.to_csv(OUTPUT_DIR / "feature_rank_stability.csv", index=False)
feature_specificity_results.to_csv(
    OUTPUT_DIR / "feature_specificity.csv", index=False
)

configuration = pd.DataFrame(
    [{
        "bootstrap_iterations": args.bootstrap_iterations,
        "rank_bootstrap_iterations": args.rank_bootstrap_iterations,
        "bootstrap_seed": args.bootstrap_seed,
        "max_missing_fraction": args.max_missing_fraction,
        "min_group_n": args.min_group_n,
        "hedges_g_threshold": args.effect_threshold,
        "kw_epsilon_squared_threshold": args.kw_effect_threshold,
        "mad_threshold": args.mad_threshold,
        "top_k": args.top_k,
        "specificity_significance_rule": "global_pairwise_q_fdr<0.05",
        "specificity_effect_rule": f"abs(hedges_g)>={args.effect_threshold}",
    }]
)
configuration.to_csv(OUTPUT_DIR / "analysis_configuration.csv", index=False)

print(f"\nSaved files to {OUTPUT_DIR}:")
print("1. stats_healthy_vs_disease.csv")
print("2. stats_each_disease_vs_healthy.csv")
print("3. stats_kruskal_all_groups.csv")
print("4. significant_healthy_vs_disease.csv")
print("5. significant_each_disease_vs_healthy.csv")
print("6. significant_kruskal_all_groups.csv")
print("7. high_priority_healthy_vs_disease.csv")
print("8. high_priority_each_disease_vs_healthy.csv")
print("9. high_priority_kruskal_all_groups.csv")
print("10. feature_reliability.csv")
print("11. feature_rank_stability.csv")
print("12. analysis_configuration.csv")
print("13. feature_specificity.csv")
