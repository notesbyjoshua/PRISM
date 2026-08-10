"""
Improved Phenotype Difference Score (PDS) Calculator
====================================================

Edit the CONFIG section, then press Run.

This version adds:
1. Stronger small-sample penalty
2. Bootstrap CI precision as part of reliability
3. Signed PDS in [-100, +100]
4. ML-ready normalized PDS in [-1, +1]

Install once:
    pip install pandas numpy scipy statsmodels
"""

# ============================================================
# CONFIG — EDIT THESE
# ============================================================

DISEASE_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_disease.csv"
HEALTHY_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_healthy.csv"
OUTPUT_FOLDER = r"pds_results_v2"
SYNDROME_COLUMN = "syndrome_name"
SINGLE_DISEASE_NAME = "Disease"
N_BOOTSTRAPS = 1000
RANDOM_STATE = 42
MIN_NONMISSING_PER_GROUP = 5

# ------------------------------------------------------------
# PDS PARAMETERS
# ------------------------------------------------------------
# Effect component:
# E = tanh(|g| / EFFECT_SCALE)

EFFECT_SCALE = 2.0

# Sample reliability:
# N = 1 - exp(-n_min / SAMPLE_SCALE)
# Larger denominator = more conservative sample-size credit

SAMPLE_SCALE = 75.0

# Reliability weights must sum to 1.0
BOOTSTRAP_DIRECTION_WEIGHT = 0.50
SAMPLE_SIZE_WEIGHT = 0.25
CI_PRECISION_WEIGHT = 0.25

# CI precision:
# C = 1 / (1 + CI_width / CI_WIDTH_SCALE)
#
# Smaller CI width -> C closer to 1
# Larger CI width -> C closer to 0

CI_WIDTH_SCALE = 1.0

# Additional hard penalty for very small disease groups
# This multiplies the final PDS.
SMALL_SAMPLE_MULTIPLIERS = {"lt_10": 0.60, "10_to_19": 0.80, "20_plus": 1.00}

# q-value multiplier
Q_MULTIPLIER_STRONG = 1.00  # q < 0.01
Q_MULTIPLIER_SIGNIFICANT = 0.80  # 0.01 <= q < 0.05
Q_MULTIPLIER_NONSIGNIFICANT = 0.25  # q >= 0.05

# ------------------------------------------------------------
# EXCLUDED NON-FACIAL COLUMNS
# ------------------------------------------------------------

EXCLUDED_COLUMNS = {
    "image_id",
    "facekit_image_id",
    "patient_id",
    "patient_name",
    "syndrome_id",
    "internal_syndrome_id",
    "syndrome_name",
    "internal_syndrome_name",
    "disease",
    "diagnosis",
    "label",
    "gene_names",
    "gene_entrez_ids",
    "hgvs",
    "hpo_terms",
    "present_features",
    "absent_features",
    "age_year",
    "age_month",
    "age_note",
    "gender",
    "ethnicity_category",
    "ethnicity_sub_category",
    "test_type",
    "OMIM",
    "train",
    "frontal_ok",
    "pose_yaw",
    "pose_pitch",
    "pose_roll",
}

# ============================================================
# IMPORTS
# ============================================================

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ============================================================
# STATISTICAL FUNCTIONS
# ============================================================


def hedges_g(disease_values, healthy_values):
    x = np.asarray(disease_values, dtype=float)
    y = np.asarray(healthy_values, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    n1 = len(x)
    n2 = len(y)
    if n1 < 2 or n2 < 2:
        return np.nan
    var1 = np.var(x, ddof=1)
    var2 = np.var(y, ddof=1)
    df = n1 + n2 - 2
    if df <= 0:
        return np.nan
    pooled_variance = ((n1 - 1) * var1 + (n2 - 1) * var2) / df
    if pooled_variance <= 0 or not np.isfinite(pooled_variance):
        return np.nan
    pooled_sd = np.sqrt(pooled_variance)
    cohens_d = (np.mean(x) - np.mean(y)) / pooled_sd
    correction = 1 - (3 / (4 * df - 1))
    return float(correction * cohens_d)


def bootstrap_g_and_stability(disease_values, healthy_values, observed_g, rng):
    """
    Returns:
        directional_stability
        bootstrap_ci_lower
        bootstrap_ci_upper
    """
    disease_values = np.asarray(disease_values, dtype=float)
    healthy_values = np.asarray(healthy_values, dtype=float)
    bootstrap_g = np.empty(N_BOOTSTRAPS, dtype=float)
    for i in range(N_BOOTSTRAPS):
        disease_sample = rng.choice(disease_values, size=len(disease_values), replace=True)
        healthy_sample = rng.choice(healthy_values, size=len(healthy_values), replace=True)
        bootstrap_g[i] = hedges_g(disease_sample, healthy_sample)
    bootstrap_g = bootstrap_g[np.isfinite(bootstrap_g)]
    if len(bootstrap_g) == 0:
        return np.nan, np.nan, np.nan
    if observed_g > 0:
        stability = np.mean(bootstrap_g > 0)
    elif observed_g < 0:
        stability = np.mean(bootstrap_g < 0)
    else:
        stability = 0.5
    ci_lower = np.percentile(bootstrap_g, 2.5)
    ci_upper = np.percentile(bootstrap_g, 97.5)
    return (float(stability), float(ci_lower), float(ci_upper))


# ============================================================
# PDS COMPONENTS
# ============================================================


def effect_component(g):
    """
    Smooth effect magnitude term.

    E = tanh(|g| / EFFECT_SCALE)
    """
    if not np.isfinite(g):
        return np.nan
    return float(np.tanh(abs(g) / EFFECT_SCALE))


def sample_reliability(n_disease, n_healthy):
    """
    Uses the smaller group as the limiting sample size.

    N = 1 - exp(-n_min / SAMPLE_SCALE)
    """
    n_min = min(n_disease, n_healthy)
    return float(1 - np.exp(-n_min / SAMPLE_SCALE))


def small_sample_multiplier(n_disease, n_healthy):
    """
    Extra safeguard against very small groups.
    """
    n_min = min(n_disease, n_healthy)
    if n_min < 10:
        return SMALL_SAMPLE_MULTIPLIERS["lt_10"]
    if n_min < 20:
        return SMALL_SAMPLE_MULTIPLIERS["10_to_19"]
    return SMALL_SAMPLE_MULTIPLIERS["20_plus"]


def ci_precision(ci_lower, ci_upper):
    """
    Precision based on bootstrap 95% CI width.

    C = 1 / (1 + CI_width / CI_WIDTH_SCALE)

    Examples if CI_WIDTH_SCALE = 1:
        width 0.25 -> C = 0.80
        width 0.50 -> C = 0.67
        width 1.00 -> C = 0.50
        width 2.00 -> C = 0.33
    """
    if not np.isfinite(ci_lower) or not np.isfinite(ci_upper):
        return np.nan
    width = max(0.0, ci_upper - ci_lower)
    return float(1 / (1 + width / CI_WIDTH_SCALE))


def q_multiplier(q):
    if not np.isfinite(q):
        return 0.0
    if q < 0.01:
        return Q_MULTIPLIER_STRONG
    if q < 0.05:
        return Q_MULTIPLIER_SIGNIFICANT
    return Q_MULTIPLIER_NONSIGNIFICANT


def pds_category(absolute_pds):
    if not np.isfinite(absolute_pds):
        return "NA"
    if absolute_pds < 20:
        return "minimal"
    if absolute_pds < 40:
        return "weak"
    if absolute_pds < 60:
        return "moderate"
    if absolute_pds < 75:
        return "strong"
    if absolute_pds < 90:
        return "very_strong"
    return "exceptional"


# ============================================================
# FEATURE DETECTION
# ============================================================


def detect_features(disease_df, healthy_df):
    features = []
    for column in disease_df.columns:
        if column not in healthy_df.columns or column in EXCLUDED_COLUMNS:
            continue
        disease_numeric = pd.to_numeric(disease_df[column], errors="coerce")
        healthy_numeric = pd.to_numeric(healthy_df[column], errors="coerce")
        if disease_numeric.notna().sum() > 0 and healthy_numeric.notna().sum() > 0:
            features.append(column)
    return features


# ============================================================
# LOAD DATA
# ============================================================

print("Loading datasets...")
disease = pd.read_csv(DISEASE_CSV)
healthy = pd.read_csv(HEALTHY_CSV)

# Keep frontal images only, if available
if "frontal_ok" in disease.columns:
    disease = disease[disease["frontal_ok"] == True].copy()
if "frontal_ok" in healthy.columns:
    healthy = healthy[healthy["frontal_ok"] == True].copy()
features = detect_features(disease, healthy)
print(f"Disease rows: {len(disease):,}")
print(f"Healthy rows: {len(healthy):,}")
print(f"Facial features: {len(features)}")

# ============================================================
# DISEASE GROUPS
# ============================================================

if SYNDROME_COLUMN is None:
    disease_groups = [(SINGLE_DISEASE_NAME, disease.copy())]
else:
    if SYNDROME_COLUMN not in disease.columns:
        raise KeyError(f"Missing column: " f"{SYNDROME_COLUMN}")
    disease[SYNDROME_COLUMN] = disease[SYNDROME_COLUMN].astype("string").str.strip()
    disease_groups = [
        (syndrome, group.copy())
        for syndrome, group in disease.groupby(SYNDROME_COLUMN, observed=True)
        if (pd.notna(syndrome) and str(syndrome).strip())
    ]
print(f"Diseases to analyze: " f"{len(disease_groups)}")

# ============================================================
# PASS 1 — RAW STATISTICS
# ============================================================

rows = []
for disease_index, (syndrome, disease_group) in enumerate(disease_groups, start=1):
    print(f"[{disease_index}/" f"{len(disease_groups)}] " f"{syndrome}")
    for feature in features:
        disease_values = pd.to_numeric(disease_group[feature], errors="coerce").dropna().to_numpy()
        healthy_values = pd.to_numeric(healthy[feature], errors="coerce").dropna().to_numpy()
        if len(disease_values) < MIN_NONMISSING_PER_GROUP or len(healthy_values) < MIN_NONMISSING_PER_GROUP:
            continue
        try:
            u_statistic, p_value = mannwhitneyu(disease_values, healthy_values, alternative="two-sided")
        except ValueError:
            u_statistic = np.nan
            p_value = 1.0
        g = hedges_g(disease_values, healthy_values)
        rows.append(
            {
                "syndrome_name": syndrome,
                "feature": feature,
                "n_disease": len(disease_values),
                "n_healthy": len(healthy_values),
                "disease_mean": float(np.mean(disease_values)),
                "healthy_mean": float(np.mean(healthy_values)),
                "disease_median": float(np.median(disease_values)),
                "healthy_median": float(np.median(healthy_values)),
                "mean_difference": float(np.mean(disease_values) - np.mean(healthy_values)),
                "hedges_g": g,
                "abs_hedges_g": (abs(g) if np.isfinite(g) else np.nan),
                "mann_whitney_u": float(u_statistic) if np.isfinite(u_statistic) else np.nan,
                "p_value": float(p_value),
            }
        )
results = pd.DataFrame(rows)
if results.empty:
    raise ValueError("No valid comparisons were produced.")

# ============================================================
# FDR CORRECTION
# ============================================================

_, global_q, _, _ = multipletests(results["p_value"], method="fdr_bh")
results["global_fdr_q_value"] = global_q
results["global_significant_fdr"] = results["global_fdr_q_value"] < 0.05
results["within_disease_fdr_q_value"] = np.nan
for syndrome in results["syndrome_name"].dropna().unique():
    mask = results["syndrome_name"] == syndrome
    _, within_q, _, _ = multipletests(results.loc[mask, "p_value"], method="fdr_bh")
    results.loc[mask, "within_disease_fdr_q_value"] = within_q

# ============================================================
# PASS 2 — BOOTSTRAP + IMPROVED PDS
# ============================================================

print("\nCalculating bootstrap " "stability and improved PDS...")
rng = np.random.default_rng(RANDOM_STATE)
extra_rows = []
for row_index, row in results.iterrows():
    if row_index == 0 or (row_index + 1) % 250 == 0 or row_index + 1 == len(results):
        print(f"PDS progress: " f"{row_index + 1:,}/" f"{len(results):,}")
    syndrome = row["syndrome_name"]
    feature = row["feature"]
    if SYNDROME_COLUMN is None:
        disease_group = disease
    else:
        disease_group = disease[disease[SYNDROME_COLUMN] == syndrome]
    disease_values = pd.to_numeric(disease_group[feature], errors="coerce").dropna().to_numpy()
    healthy_values = pd.to_numeric(healthy[feature], errors="coerce").dropna().to_numpy()
    g = row["hedges_g"]
    stability, ci_lower, ci_upper = bootstrap_g_and_stability(disease_values, healthy_values, g, rng)
    E = effect_component(g)
    N = sample_reliability(len(disease_values), len(healthy_values))
    C = ci_precision(ci_lower, ci_upper)
    R = BOOTSTRAP_DIRECTION_WEIGHT * stability + SAMPLE_SIZE_WEIGHT * N + CI_PRECISION_WEIGHT * C
    Q = q_multiplier(row["global_fdr_q_value"])
    small_n_multiplier = small_sample_multiplier(len(disease_values), len(healthy_values))
    pds_normalized = np.sign(g) * E * R * Q * small_n_multiplier
    pds_signed = 100 * pds_normalized
    extra_rows.append(
        {
            "bootstrap_direction_stability": stability,
            "hedges_g_bootstrap_ci_lower": ci_lower,
            "hedges_g_bootstrap_ci_upper": ci_upper,
            "hedges_g_bootstrap_ci_width": (
                ci_upper - ci_lower if (np.isfinite(ci_upper) and np.isfinite(ci_lower)) else np.nan
            ),
            "pds_effect_component": E,
            "pds_sample_reliability": N,
            "pds_ci_precision": C,
            "pds_reliability_component": R,
            "pds_q_multiplier": Q,
            "pds_small_sample_multiplier": small_n_multiplier,
            "pds_normalized": float(pds_normalized),
            "pds_signed": float(pds_signed),
            "pds_absolute": float(abs(pds_signed)),
        }
    )
results = pd.concat([results.reset_index(drop=True), pd.DataFrame(extra_rows)], axis=1)

# ============================================================
# PDS LABELS + PERCENTILES
# ============================================================

results["pds_category"] = results["pds_absolute"].apply(pds_category)
results["pds_within_disease_percentile"] = (
    results.groupby("syndrome_name", observed=True)["pds_absolute"].rank(pct=True, method="average") * 100
)
results = results.sort_values(["syndrome_name", "pds_absolute"], ascending=[True, False])

# ============================================================
# SAVE OUTPUTS
# ============================================================

output = Path(OUTPUT_FOLDER)
output.mkdir(parents=True, exist_ok=True)
results.to_csv(output / "pds_v2_all_disease_feature_results.csv", index=False)
results[results["global_significant_fdr"]].to_csv(output / "pds_v2_significant_features.csv", index=False)
(results.groupby("syndrome_name", observed=True).head(15).to_csv(output / "pds_v2_top_15_per_disease.csv", index=False))
results.pivot(index="syndrome_name", columns="feature", values="pds_signed").to_csv(output / "pds_v2_signed_matrix.csv")
results.pivot(index="syndrome_name", columns="feature", values="pds_normalized").to_csv(
    output / "pds_v2_normalized_matrix_for_ml.csv"
)

# ============================================================
# SAVE CONFIGURATION
# ============================================================

configuration = {
    "formula": (
        "sign(g) * "
        "tanh(abs(g)/2) * "
        "(0.50*bootstrap_direction_stability + "
        "0.25*sample_reliability + "
        "0.25*ci_precision) * "
        "q_multiplier * "
        "small_sample_multiplier"
    ),
    "sample_reliability": ("1 - exp(-min(n_disease,n_healthy)/75)"),
    "ci_precision": ("1 / (1 + CI_width / 1.0)"),
    "q_multiplier": {"q < 0.01": 1.0, "0.01 <= q < 0.05": 0.8, "q >= 0.05": 0.25},
    "small_sample_multiplier": {"n < 10": 0.60, "10 <= n < 20": 0.80, "n >= 20": 1.00},
    "N_BOOTSTRAPS": N_BOOTSTRAPS,
    "EFFECT_SCALE": EFFECT_SCALE,
    "SAMPLE_SCALE": SAMPLE_SCALE,
    "CI_WIDTH_SCALE": CI_WIDTH_SCALE,
}
with open(output / "pds_v2_configuration.json", "w", encoding="utf-8") as file:
    json.dump(configuration, file, indent=2)

# ============================================================
# FINISH
# ============================================================

print("\nDone.")
print(f"Comparisons: " f"{len(results):,}")
print(f"Max absolute PDS: " f"{results['pds_absolute'].max():.2f}")
print(f"Median absolute PDS: " f"{results['pds_absolute'].median():.2f}")
print(f"PDS >= 60: " f"{(results['pds_absolute'] >= 60).sum():,}")
print(f"PDS >= 75: " f"{(results['pds_absolute'] >= 75).sum():,}")
print(f"PDS >= 90: " f"{(results['pds_absolute'] >= 90).sum():,}")
print("\nSaved to:")
print(output.resolve())
