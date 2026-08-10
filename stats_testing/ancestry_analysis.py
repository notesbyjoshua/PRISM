"""
Two-Layer Race/Ethnicity Analysis for FaceKit
==============================================

This script harmonizes GMDB disease metadata and FairFace healthy metadata
into TWO analysis layers, then runs Analysis A and Analysis B for each layer.

LAYER 1 — BROAD ETHNICITY
Proper labels:
    European
    African
    Asian
    Others

GMDB:
    European -> European
    African   -> African
    Asian     -> Asian
    Others    -> Others
    Unknown   -> excluded

FairFace:
    Black -> African
    White -> European
    East Asian / Central Asian / West Asian / Middle Eastern /
    Indian / Southeast Asian -> Asian
    Latino_Hispanic / Latino-Hispanic / Others -> Others
    Unknown -> excluded


LAYER 2 — SUBCATEGORY
Proper labels:
    African
    East Asian
    South Asian
    SE Asian
    Latino/Hispanic
    Middle Eastern
    European

GMDB:
    ethnicity_category == African -> African
        (subcategory ignored)

    ethnicity_category == European -> European
        (subcategory ignored)

    ethnicity_sub_category:
        East Asian -> East Asian
        South Asian -> South Asian
        SE Asian / Southeast Asian -> SE Asian
        American - Latin/Hispanic -> Latino/Hispanic
        Middle-East / Middle Eastern / West Asian -> Middle Eastern

    Asian Others -> excluded
    ethnicity_category == Others without a usable mapped subcategory -> excluded
    Unknown -> excluded

FairFace:
    Black -> African
    White -> European
    East Asian -> East Asian
    Southeast Asian / SE Asian -> SE Asian
    Indian -> South Asian
    Middle Eastern -> Middle Eastern
    Latino_Hispanic / Latino-Hispanic -> Latino/Hispanic
    Unknown / Others -> excluded


ANALYSIS A
----------
Healthy controls only:
    1. Kruskal-Wallis across all eligible groups
    2. Pairwise Mann-Whitney U between groups
    3. Hedges' g
    4. BH-FDR correction

Question:
    Which FaceKit measurements normally differ across race/ethnicity groups?


ANALYSIS B
----------
For each syndrome × race/ethnicity group:
    disease patients in group X
        vs
    healthy FairFace controls in matched group X

Calculates:
    Mann-Whitney U
    Hedges' g
    bootstrap 95% CI
    within-group BH-FDR
    global BH-FDR
    high-priority flag

Question:
    Within a matched race/ethnicity group, which facial measurements
    differ between a syndrome and healthy controls?


IMPORTANT TERMINOLOGY
---------------------
These are harmonized dataset race/ethnicity labels, NOT genetically inferred
ancestry. In a paper/poster, describe them as:
    "harmonized race/ethnicity groups"

Edit CONFIG and press Run.

Install once:
    pip install pandas numpy scipy statsmodels
"""

# ============================================================
# CONFIG — EDIT THESE
# ============================================================

FAIRFACE_HEALTHY_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_healthy.csv"
GMDB_DISEASE_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_disease.csv"
OUTPUT_FOLDER = r"stats_results/ethnicity_analysis"

FAIRFACE_RACE_COLUMN = "race"
GMDB_ETHNICITY_COLUMN = "image_dataset_ethnicity_category"
GMDB_SUBCATEGORY_COLUMN = "image_dataset_ethnicity_sub_category"
SYNDROME_COLUMN = "internal_syndrome_name"

# Exact FaceKit feature range
FIRST_FACEKIT_FEATURE = "eb_thickness_r"
LAST_FACEKIT_FEATURE = "cheek_area_asym"

# Sample-size thresholds
MIN_HEALTHY_PER_GROUP = 20
MIN_DISEASE_PER_GROUP = 2
MIN_MATCHED_HEALTHY = 20

# Multiple testing
FDR_ALPHA = 0.05

# Bootstrap
# 100 for quick testing; 1000 for final analysis
N_BOOTSTRAPS = 1000
RANDOM_STATE = 42

# Require frontal_ok when the column exists
REQUIRE_FRONTAL = True

# ============================================================
# IMPORTS
# ============================================================

import itertools
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ============================================================
# STRING NORMALIZATION
# ============================================================


def norm_label(value):
    """
    Normalize labels so differences in capitalization, underscores,
    hyphens, and repeated spaces do not break mapping.
    """
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    text = text.lower()

    # Make separators comparable
    text = text.replace("_", " ")
    text = text.replace("/", " ")
    text = text.replace("-", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if text in {"unknown", "unk", "nan", "none", "not known", "unspecified", "n a", "na"}:
        return None
    return text


# ============================================================
# LAYER 1 MAPPING
# ============================================================


def map_gmdb_layer1(value):
    """
    GMDB broad ethnicity -> Layer 1.
    """
    x = norm_label(value)
    if x is None:
        return np.nan
    if x in {"european", "europe", "caucasian", "white"}:
        return "European"
    if x in {"african", "africa", "black"}:
        return "African"
    if x in {"asian", "asia"}:
        return "Asian"
    if x in {"others", "other"}:
        return "Others"
    return np.nan


def map_fairface_layer1(value):
    """
    FairFace race -> Layer 1.
    """
    x = norm_label(value)
    if x is None:
        return np.nan
    if x in {"black", "african"}:
        return "African"
    if x in {"white", "european"}:
        return "European"

    # User-specified Asian umbrella
    if x in {
        "east asian",
        "central asian",
        "west asian",
        "western asian",
        "middle eastern",
        "middle east",
        "indian",
        "south asian",
        "southeast asian",
        "south east asian",
        "se asian",
        "asian",
    }:
        return "Asian"
    if x in {"latino hispanic", "latino", "hispanic", "other", "others"}:
        return "Others"
    return np.nan


# ============================================================
# LAYER 2 MAPPING
# ============================================================


def map_gmdb_layer2(ethnicity, subcategory):
    """
    GMDB category + subcategory -> Layer 2.

    Priority:
        African broad category -> African
        European broad category -> European
        Otherwise use valid subcategory mapping.
    """
    broad = norm_label(ethnicity)
    sub = norm_label(subcategory)
    if broad is None:
        return np.nan

    # African: do not differentiate subcategory
    if broad in {"african", "africa", "black"}:
        return "African"

    # European has no useful subcategory in GMDB
    if broad in {"european", "europe", "caucasian", "white"}:
        return "European"

    # Unknown broad ethnicity already removed above.
    # "Others" is allowed only if its subcategory maps to a valid target.
    if sub is None:
        return np.nan
    if sub in {"east asian", "eastern asian"}:
        return "East Asian"
    if sub in {"south asian", "southern asian", "indian"}:
        return "South Asian"
    if sub in {"se asian", "southeast asian", "south east asian", "south eastern asian"}:
        return "SE Asian"
    if sub in {
        "american latin hispanic",
        "american latino hispanic",
        "latin hispanic",
        "latino hispanic",
        "latino",
        "hispanic",
    }:
        return "Latino/Hispanic"
    if sub in {"middle east", "middle eastern", "middle east west asian", "west asian", "western asian"}:
        return "Middle Eastern"

    # Explicitly excluded
    if sub in {"asian others", "asian other", "other asian", "others asian"}:
        return np.nan
    return np.nan


def map_fairface_layer2(value):
    """
    FairFace race -> Layer 2.
    """
    x = norm_label(value)
    if x is None:
        return np.nan
    if x in {"black", "african"}:
        return "African"
    if x in {"white", "european"}:
        return "European"
    if x in {"east asian", "eastern asian"}:
        return "East Asian"
    if x in {"southeast asian", "south east asian", "se asian"}:
        return "SE Asian"
    if x in {"indian", "south asian", "southern asian"}:
        return "South Asian"
    if x in {"middle eastern", "middle east", "west asian", "western asian"}:
        return "Middle Eastern"
    if x in {"latino hispanic", "latino", "hispanic"}:
        return "Latino/Hispanic"

    # User requested "others" excluded from Layer 2.
    if x in {"other", "others", "asian others", "asian other"}:
        return np.nan
    return np.nan


# ============================================================
# STATISTICS
# ============================================================


def hedges_g(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    n1 = len(x)
    n2 = len(y)
    if n1 < 2 or n2 < 2:
        return np.nan
    v1 = np.var(x, ddof=1)
    v2 = np.var(y, ddof=1)
    df = n1 + n2 - 2
    if df <= 0:
        return np.nan
    pooled_variance = (((n1 - 1) * v1) + ((n2 - 1) * v2)) / df
    if pooled_variance <= 0 or not np.isfinite(pooled_variance):
        return np.nan
    d = (np.mean(x) - np.mean(y)) / np.sqrt(pooled_variance)
    correction = 1 - 3 / (4 * df - 1)
    return float(correction * d)


def effect_category(g):
    if not np.isfinite(g):
        return "NA"
    a = abs(g)
    if a < 0.2:
        return "negligible"
    if a < 0.5:
        return "small"
    if a < 0.8:
        return "medium"
    return "large"


def bootstrap_hedges_g_ci(disease_values, healthy_values, rng):
    disease_values = np.asarray(disease_values, dtype=float)
    healthy_values = np.asarray(healthy_values, dtype=float)
    estimates = np.empty(N_BOOTSTRAPS, dtype=float)
    for i in range(N_BOOTSTRAPS):
        d_boot = rng.choice(disease_values, size=len(disease_values), replace=True)
        h_boot = rng.choice(healthy_values, size=len(healthy_values), replace=True)
        estimates[i] = hedges_g(d_boot, h_boot)
    estimates = estimates[np.isfinite(estimates)]
    if len(estimates) == 0:
        return (np.nan, np.nan)
    return (float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5)))


def add_global_fdr(dataframe, p_column="p_value"):
    dataframe = dataframe.copy()
    if len(dataframe) == 0:
        return dataframe
    valid = dataframe[p_column].notna()
    dataframe["global_fdr_q_value"] = np.nan
    dataframe["global_significant_fdr"] = False
    if valid.sum() > 0:
        reject, q, _, _ = multipletests(dataframe.loc[valid, p_column], alpha=FDR_ALPHA, method="fdr_bh")
        dataframe.loc[valid, "global_fdr_q_value"] = q
        dataframe.loc[valid, "global_significant_fdr"] = reject
    return dataframe


# ============================================================
# LOAD DATA
# ============================================================

print("Loading FairFace healthy data...")
healthy = pd.read_csv(FAIRFACE_HEALTHY_CSV)
print("Loading GMDB disease data...")
disease = pd.read_csv(GMDB_DISEASE_CSV)

# ============================================================
# VALIDATE REQUIRED COLUMNS
# ============================================================

required_healthy = {FAIRFACE_RACE_COLUMN, FIRST_FACEKIT_FEATURE, LAST_FACEKIT_FEATURE}
required_disease = {
    GMDB_ETHNICITY_COLUMN,
    GMDB_SUBCATEGORY_COLUMN,
    SYNDROME_COLUMN,
    FIRST_FACEKIT_FEATURE,
    LAST_FACEKIT_FEATURE,
}
missing_healthy = required_healthy - set(healthy.columns)
missing_disease = required_disease - set(disease.columns)
if missing_healthy:
    raise KeyError("Healthy CSV missing: " f"{sorted(missing_healthy)}")
if missing_disease:
    raise KeyError("Disease CSV missing: " f"{sorted(missing_disease)}")

# ============================================================
# FRONTAL FILTER
# ============================================================

if REQUIRE_FRONTAL:
    if "frontal_ok" in healthy.columns:
        healthy = healthy[healthy["frontal_ok"] == True].copy()
    if "frontal_ok" in disease.columns:
        disease = disease[disease["frontal_ok"] == True].copy()

# ============================================================
# CLEAN SYNDROME LABEL
# ============================================================

disease[SYNDROME_COLUMN] = disease[SYNDROME_COLUMN].astype("string").str.strip()
disease = disease[disease[SYNDROME_COLUMN].notna()].copy()

# ============================================================
# CREATE HARMONIZED LAYERS
# ============================================================

print("\nCreating harmonized Layer 1...")
healthy["layer1_group"] = healthy[FAIRFACE_RACE_COLUMN].apply(map_fairface_layer1)
disease["layer1_group"] = disease[GMDB_ETHNICITY_COLUMN].apply(map_gmdb_layer1)
print("Creating harmonized Layer 2...")
healthy["layer2_group"] = healthy[FAIRFACE_RACE_COLUMN].apply(map_fairface_layer2)
disease["layer2_group"] = [
    map_gmdb_layer2(ethnicity, subcategory)
    for ethnicity, subcategory in zip(disease[GMDB_ETHNICITY_COLUMN], disease[GMDB_SUBCATEGORY_COLUMN])
]

# ============================================================
# EXACT 125 FACEKIT FEATURES
# ============================================================

healthy_features = list(healthy.loc[:, FIRST_FACEKIT_FEATURE:LAST_FACEKIT_FEATURE].columns)
disease_features = list(disease.loc[:, FIRST_FACEKIT_FEATURE:LAST_FACEKIT_FEATURE].columns)
feature_columns = [feature for feature in healthy_features if feature in disease_features]
if len(feature_columns) != 125:
    raise ValueError(f"Expected exactly 125 shared FaceKit features, " f"found {len(feature_columns)}.")
for feature in feature_columns:
    healthy[feature] = pd.to_numeric(healthy[feature], errors="coerce")
    disease[feature] = pd.to_numeric(disease[feature], errors="coerce")
print("FaceKit features:", len(feature_columns))

# ============================================================
# OUTPUT
# ============================================================

output = Path(OUTPUT_FOLDER)
output.mkdir(parents=True, exist_ok=True)

# ============================================================
# SAVE HARMONIZATION AUDIT TABLES
# ============================================================
# FairFace original -> Layer 1 / Layer 2

fairface_mapping_audit = (
    healthy[[FAIRFACE_RACE_COLUMN, "layer1_group", "layer2_group"]].value_counts(dropna=False).reset_index(name="n")
)
fairface_mapping_audit.to_csv(output / "fairface_harmonization_audit.csv", index=False)

# GMDB original -> Layer 1 / Layer 2
gmdb_mapping_audit = (
    disease[[GMDB_ETHNICITY_COLUMN, GMDB_SUBCATEGORY_COLUMN, "layer1_group", "layer2_group"]]
    .value_counts(dropna=False)
    .reset_index(name="n")
)
gmdb_mapping_audit.to_csv(output / "gmdb_harmonization_audit.csv", index=False)


# ============================================================
# GENERIC ANALYSIS FUNCTION
# ============================================================


def run_layer_analysis(layer_name, group_column):
    print("\n========================================")
    print(f"{layer_name.upper()}")
    print("========================================")
    layer_folder = output / layer_name
    layer_folder.mkdir(parents=True, exist_ok=True)
    healthy_layer = healthy[healthy[group_column].notna()].copy()
    disease_layer = disease[disease[group_column].notna()].copy()

    # --------------------------------------------------------
    # COUNTS
    # --------------------------------------------------------

    healthy_counts = healthy_layer[group_column].value_counts().rename_axis("group").reset_index(name="healthy_images")
    healthy_counts.to_csv(layer_folder / "healthy_group_counts.csv", index=False)
    syndrome_counts = (
        disease_layer.groupby([SYNDROME_COLUMN, group_column], observed=True)
        .size()
        .reset_index(name="disease_images")
        .rename(columns={group_column: "group"})
    )
    syndrome_counts.to_csv(layer_folder / "syndrome_group_counts.csv", index=False)
    print("Healthy counts:")
    print(healthy_counts.to_string(index=False))

    # ========================================================
    # ANALYSIS A
    # ========================================================

    print("\nAnalysis A: healthy baseline")
    healthy_group_sizes = healthy_layer[group_column].value_counts()
    eligible_groups = healthy_group_sizes[healthy_group_sizes >= MIN_HEALTHY_PER_GROUP].index.tolist()
    print("Eligible healthy groups:", eligible_groups)
    if len(eligible_groups) < 2:
        print("WARNING: fewer than 2 eligible healthy groups. " "Skipping Analysis A.")
        kw = pd.DataFrame()
        pairwise = pd.DataFrame()
    else:

        # ----------------------------------------------------
        # A1. KRUSKAL-WALLIS
        # ----------------------------------------------------

        kw_rows = []
        for feature_index, feature in enumerate(feature_columns, start=1):
            groups = []
            group_ns = {}
            group_means = {}
            group_medians = {}
            for group_name in eligible_groups:
                values = healthy_layer.loc[healthy_layer[group_column] == group_name, feature].dropna().to_numpy()
                if len(values) < MIN_HEALTHY_PER_GROUP:
                    continue
                groups.append(values)
                group_ns[group_name] = len(values)
                group_means[group_name] = float(np.mean(values))
                group_medians[group_name] = float(np.median(values))
            if len(groups) < 2:
                continue
            try:
                H, p = kruskal(*groups)
            except ValueError:
                H = 0.0
                p = 1.0
            kw_rows.append(
                {
                    "feature": feature,
                    "kruskal_H": float(H),
                    "p_value": float(p),
                    "number_of_groups": len(groups),
                    "group_sample_sizes": "|".join(f"{k}:{v}" for k, v in group_ns.items()),
                    "group_means": "|".join(f"{k}:{v:.6g}" for k, v in group_means.items()),
                    "group_medians": "|".join(f"{k}:{v:.6g}" for k, v in group_medians.items()),
                }
            )
            if feature_index == 1 or feature_index % 25 == 0 or feature_index == len(feature_columns):
                print(f"  Kruskal-Wallis: " f"{feature_index}/" f"{len(feature_columns)}")
        kw = pd.DataFrame(kw_rows)
        if len(kw) > 0:
            reject, q, _, _ = multipletests(kw["p_value"], alpha=FDR_ALPHA, method="fdr_bh")
            kw["fdr_q_value"] = q
            kw["significant_fdr"] = reject
            kw.to_csv(layer_folder / "analysis_A_healthy_kruskal_all_groups.csv", index=False)
            kw[kw["significant_fdr"]].to_csv(layer_folder / "analysis_A_healthy_kruskal_significant.csv", index=False)

        # ----------------------------------------------------
        # A2. PAIRWISE HEALTHY COMPARISONS
        # ----------------------------------------------------

        pair_rows = []
        group_pairs = list(itertools.combinations(eligible_groups, 2))
        for pair_index, (group_1, group_2) in enumerate(group_pairs, start=1):
            print(f"  Pair " f"{pair_index}/" f"{len(group_pairs)}: " f"{group_1} vs " f"{group_2}")
            for feature in feature_columns:
                x = healthy_layer.loc[healthy_layer[group_column] == group_1, feature].dropna().to_numpy()
                y = healthy_layer.loc[healthy_layer[group_column] == group_2, feature].dropna().to_numpy()
                if len(x) < MIN_HEALTHY_PER_GROUP or len(y) < MIN_HEALTHY_PER_GROUP:
                    continue
                try:
                    U, p = mannwhitneyu(x, y, alternative="two-sided")
                except ValueError:
                    U = np.nan
                    p = 1.0
                g = hedges_g(x, y)
                pair_rows.append(
                    {
                        "group_1": group_1,
                        "group_2": group_2,
                        "feature": feature,
                        "n_group_1": len(x),
                        "n_group_2": len(y),
                        "mean_group_1": float(np.mean(x)),
                        "mean_group_2": float(np.mean(y)),
                        "median_group_1": float(np.median(x)),
                        "median_group_2": float(np.median(y)),
                        "mann_whitney_u": float(U) if np.isfinite(U) else np.nan,
                        "p_value": float(p),
                        "hedges_g": g,
                        "abs_hedges_g": abs(g) if np.isfinite(g) else np.nan,
                        "effect_category": effect_category(g),
                        "direction": (
                            f"higher_in_{group_1}"
                            if (np.isfinite(g) and g > 0)
                            else f"higher_in_{group_2}" if (np.isfinite(g) and g < 0) else "no_direction"
                        ),
                    }
                )
        pairwise = pd.DataFrame(pair_rows)
        if len(pairwise) > 0:
            pairwise["within_pair_fdr_q_value"] = np.nan
            pairwise["within_pair_significant_fdr"] = False
            for _, indices in pairwise.groupby(["group_1", "group_2"], observed=True).groups.items():
                indices = list(indices)
                reject, q, _, _ = multipletests(pairwise.loc[indices, "p_value"], alpha=FDR_ALPHA, method="fdr_bh")
                pairwise.loc[indices, "within_pair_fdr_q_value"] = q
                pairwise.loc[indices, "within_pair_significant_fdr"] = reject
            pairwise.to_csv(layer_folder / "analysis_A_healthy_pairwise_groups.csv", index=False)
            pairwise[pairwise["within_pair_significant_fdr"]].to_csv(
                layer_folder / "analysis_A_healthy_pairwise_significant.csv", index=False
            )

    # ========================================================
    # ANALYSIS B
    # ========================================================

    print("\nAnalysis B: matched disease vs healthy")
    disease_group_counts = disease_layer.groupby([SYNDROME_COLUMN, group_column], observed=True).size()
    eligible_disease_groups = []
    for (syndrome, group_name), disease_n in disease_group_counts.items():
        healthy_n = int((healthy_layer[group_column] == group_name).sum())
        if disease_n >= MIN_DISEASE_PER_GROUP and healthy_n >= MIN_MATCHED_HEALTHY:
            eligible_disease_groups.append((syndrome, group_name, int(disease_n), healthy_n))
    print("Eligible syndrome × group comparisons:", len(eligible_disease_groups))
    rng = np.random.default_rng(RANDOM_STATE)
    b_rows = []
    for comparison_index, (syndrome, group_name, disease_n_total, healthy_n_total) in enumerate(
        eligible_disease_groups, start=1
    ):
        print(
            f"  [{comparison_index}/"
            f"{len(eligible_disease_groups)}] "
            f"{syndrome} | "
            f"{group_name} | "
            f"disease n="
            f"{disease_n_total} | "
            f"healthy n="
            f"{healthy_n_total}"
        )
        disease_group = disease_layer[
            (disease_layer[SYNDROME_COLUMN] == syndrome) & (disease_layer[group_column] == group_name)
        ]
        healthy_group = healthy_layer[healthy_layer[group_column] == group_name]
        for feature in feature_columns:
            d = disease_group[feature].dropna().to_numpy()
            h = healthy_group[feature].dropna().to_numpy()
            n_disease = len(d)
            if n_disease >= 20:
                evidence_level = "primary"
            elif n_disease >= 10:
                evidence_level = "moderate"
            elif n_disease >= 5:
                evidence_level = "exploratory"
            else:
                evidence_level = "insufficient"
            if len(d) < MIN_DISEASE_PER_GROUP or len(h) < MIN_MATCHED_HEALTHY:
                continue
            try:
                U, p = mannwhitneyu(d, h, alternative="two-sided")
            except ValueError:
                U = np.nan
                p = 1.0
            g = hedges_g(d, h)
            ci_low, ci_high = bootstrap_hedges_g_ci(d, h, rng)
            disease_mean = float(np.mean(d))
            healthy_mean = float(np.mean(h))
            disease_median = float(np.median(d))
            healthy_median = float(np.median(h))
            b_rows.append(
                {
                    "syndrome_name": syndrome,
                    "group": group_name,
                    "feature": feature,
                    "n_disease": len(d),
                    "n_healthy_matched": len(h),
                    "evidence_level": evidence_level,
                    "disease_mean": disease_mean,
                    "healthy_mean": healthy_mean,
                    "mean_difference": disease_mean - healthy_mean,
                    "disease_median": disease_median,
                    "healthy_median": healthy_median,
                    "median_difference": disease_median - healthy_median,
                    "mann_whitney_u": float(U) if np.isfinite(U) else np.nan,
                    "p_value": float(p),
                    "hedges_g": g,
                    "abs_hedges_g": abs(g) if np.isfinite(g) else np.nan,
                    "hedges_g_ci_lower": ci_low,
                    "hedges_g_ci_upper": ci_high,
                    "ci_excludes_zero": bool(
                        np.isfinite(ci_low) and np.isfinite(ci_high) and (ci_low > 0 or ci_high < 0)
                    ),
                    "effect_category": effect_category(g),
                    "direction": (
                        "higher_in_disease"
                        if (np.isfinite(g) and g > 0)
                        else "lower_in_disease" if (np.isfinite(g) and g < 0) else "no_direction"
                    ),
                }
            )
    analysis_b = pd.DataFrame(b_rows)
    if len(analysis_b) == 0:
        print("WARNING: Analysis B produced no comparisons.")
        return

    # --------------------------------------------------------
    # WITHIN SYNDROME × GROUP FDR
    # --------------------------------------------------------

    analysis_b["within_group_fdr_q_value"] = np.nan
    analysis_b["within_group_significant_fdr"] = False
    for _, indices in analysis_b.groupby(["syndrome_name", "group"], observed=True).groups.items():
        indices = list(indices)
        reject, q, _, _ = multipletests(analysis_b.loc[indices, "p_value"], alpha=FDR_ALPHA, method="fdr_bh")
        analysis_b.loc[indices, "within_group_fdr_q_value"] = q
        analysis_b.loc[indices, "within_group_significant_fdr"] = reject

    # --------------------------------------------------------
    # GLOBAL FDR
    # --------------------------------------------------------

    analysis_b = add_global_fdr(analysis_b, p_column="p_value")

    # --------------------------------------------------------
    # HIGH PRIORITY
    # --------------------------------------------------------

    analysis_b["high_priority"] = (
        analysis_b["within_group_significant_fdr"]
        & (analysis_b["abs_hedges_g"] >= 0.8)
        & analysis_b["ci_excludes_zero"]
        & analysis_b["evidence_level"].isin(["primary", "moderate"])
    )
    analysis_b = analysis_b.sort_values(
        ["syndrome_name", "group", "high_priority", "abs_hedges_g"], ascending=[True, True, False, False]
    )
    analysis_b.to_csv(layer_folder / "analysis_B_matched_disease_vs_healthy.csv", index=False)
    analysis_b[analysis_b["within_group_significant_fdr"]].to_csv(
        layer_folder / "analysis_B_matched_significant.csv", index=False
    )
    analysis_b[analysis_b["high_priority"]].to_csv(layer_folder / "analysis_B_matched_high_priority.csv", index=False)

    # --------------------------------------------------------
    # TOP 15
    # --------------------------------------------------------

    top15 = (
        analysis_b[analysis_b["within_group_significant_fdr"]]
        .sort_values(["syndrome_name", "group", "abs_hedges_g"], ascending=[True, True, False])
        .groupby(["syndrome_name", "group"], observed=True)
        .head(15)
    )
    top15.to_csv(layer_folder / "analysis_B_top_15_per_syndrome_group.csv", index=False)

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    summary = (
        analysis_b.groupby(["syndrome_name", "group"], observed=True)
        .agg(
            features_tested=("feature", "count"),
            significant_features=("within_group_significant_fdr", "sum"),
            high_priority_features=("high_priority", "sum"),
            max_abs_hedges_g=("abs_hedges_g", "max"),
            median_abs_hedges_g=("abs_hedges_g", "median"),
            disease_n=("n_disease", "max"),
            matched_healthy_n=("n_healthy_matched", "max"),
        )
        .reset_index()
    )
    summary.to_csv(layer_folder / "analysis_B_syndrome_group_summary.csv", index=False)
    print(f"\n{layer_name} complete:")
    if len(kw) > 0:
        print("  Analysis A significant KW features:", int(kw["significant_fdr"].sum()))
    print("  Analysis B comparisons:", f"{len(analysis_b):,}")
    print("  Analysis B significant:", f"{int(analysis_b['within_group_significant_fdr'].sum()):,}")
    print("  Analysis B high priority:", f"{int(analysis_b['high_priority'].sum()):,}")


# ============================================================
# RUN BOTH LAYERS
# ============================================================

run_layer_analysis(layer_name="layer1_broad_ethnicity", group_column="layer1_group")
run_layer_analysis(layer_name="layer2_subcategory", group_column="layer2_group")

# ============================================================
# FINAL MAPPING COUNTS
# ============================================================

mapping_summary_rows = []
for dataset_name, dataframe in [("FairFace healthy", healthy), ("GMDB disease", disease)]:
    for layer_column in ["layer1_group", "layer2_group"]:
        total = len(dataframe)
        kept = int(dataframe[layer_column].notna().sum())
        excluded = total - kept
        mapping_summary_rows.append(
            {
                "dataset": dataset_name,
                "layer": layer_column,
                "total_rows": total,
                "kept_rows": kept,
                "excluded_rows": excluded,
                "percent_kept": 100 * kept / total if total else np.nan,
            }
        )
pd.DataFrame(mapping_summary_rows).to_csv(output / "harmonization_summary.csv", index=False)
print("\n========================================")
print("ALL ANALYSES COMPLETE")
print("========================================")
print("Outputs saved to:")
print(output.resolve())
