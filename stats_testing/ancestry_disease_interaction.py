"""
Disease × Race/Ethnicity Interaction
==================================================

This script uses the SAME two-layer harmonization scheme as Ancestry Analysis.

For each:
    syndrome × FaceKit feature × pair of harmonized groups

it creates four cells:

                    Group 1             Group 2
Healthy             H1                  H2
Disease             D1                  D2

Then it fits:

    feature ~ disease + group + disease:group

The disease:group coefficient is the key Analysis C result.

Interpretation:
    interaction_beta > 0:
        the disease-associated change is more positive in Group 2 than Group 1

    interaction_beta < 0:
        the disease-associated change is more negative in Group 2 than Group 1

Mathematically:

    interaction_beta =
        (mean_D2 - mean_H2) - (mean_D1 - mean_H1)

This is a difference-in-differences interaction.

The script also reports:
    - Hedges' g within Group 1
    - Hedges' g within Group 2
    - delta_g = g_group2 - g_group1
    - robust HC3 interaction p-value
    - BH-FDR q-values
    - direction reversal flags
    - effect-concordance flags
    - evidence level based on disease sample size

IMPORTANT:
These are harmonized dataset race/ethnicity groups, not genetically inferred ancestry.

Install once:
    pip install pandas numpy scipy statsmodels
"""

# ============================================================
# CONFIG
# ============================================================

FAIRFACE_HEALTHY_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_healthy.csv"
GMDB_DISEASE_CSV = r"/Users/joshua/Documents/PRISM/data/master_dataset_disease.csv"
OUTPUT_FOLDER = r"stats_results/ethnicity_analysis"

FAIRFACE_RACE_COLUMN = "race"
GMDB_ETHNICITY_COLUMN = "image_dataset_ethnicity_category"
GMDB_SUBCATEGORY_COLUMN = "image_dataset_ethnicity_sub_category"
SYNDROME_COLUMN = "internal_syndrome_name"

FIRST_FACEKIT_FEATURE = "eb_thickness_r"
LAST_FACEKIT_FEATURE = "cheek_area_asym"

# Analysis C needs disease observations in BOTH groups.
# Keep this low enough to explore sparse groups, but use evidence_level
# when deciding what is strong enough to report.

MIN_DISEASE_PER_GROUP = 2
MIN_HEALTHY_PER_GROUP = 20

FDR_ALPHA = 0.05

# Interaction candidates can be prioritized by the difference in
# ancestry-matched Hedges' g values.
MIN_ABS_DELTA_G_FOR_PRIORITY = 0.8

REQUIRE_FRONTAL = True


# ============================================================
# IMPORTS
# ============================================================

import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# ============================================================
# LABEL HARMONIZATION
# ============================================================


def norm_label(value):
    if pd.isna(value):
        return None

    text = str(value).strip().lower()
    if not text:
        return None

    text = text.replace("_", " ").replace("/", " ").replace("-", " ")
    text = re.sub(r"\s+", " ", text).strip()

    if text in {"unknown", "unk", "nan", "none", "not known", "unspecified", "n a", "na"}:
        return None

    return text


def map_gmdb_layer1(value):
    x = norm_label(value)

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
    x = norm_label(value)

    if x in {"black", "african"}:
        return "African"
    if x in {"white", "european"}:
        return "European"
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


def map_gmdb_layer2(ethnicity, subcategory):
    broad = norm_label(ethnicity)
    sub = norm_label(subcategory)

    if broad in {"african", "africa", "black"}:
        return "African"

    if broad in {"european", "europe", "caucasian", "white"}:
        return "European"

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

    return np.nan


def map_fairface_layer2(value):
    x = norm_label(value)

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

    return np.nan


# ============================================================
# STATISTICS HELPERS
# ============================================================


def hedges_g(disease_values, healthy_values):
    x = np.asarray(disease_values, dtype=float)
    y = np.asarray(healthy_values, dtype=float)

    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]

    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return np.nan

    v1 = np.var(x, ddof=1)
    v2 = np.var(y, ddof=1)
    df = n1 + n2 - 2

    pooled_var = (((n1 - 1) * v1) + ((n2 - 1) * v2)) / df
    if pooled_var <= 0 or not np.isfinite(pooled_var):
        return np.nan

    d = (np.mean(x) - np.mean(y)) / np.sqrt(pooled_var)
    correction = 1 - 3 / (4 * df - 1)

    return float(correction * d)


def evidence_level(n1, n2):
    n_min = min(n1, n2)

    if n_min >= 20:
        return "primary"
    if n_min >= 10:
        return "moderate"
    if n_min >= 5:
        return "exploratory"

    return "insufficient"


def add_fdr_by_syndrome(data):
    data = data.copy()
    data["within_syndrome_fdr_q_value"] = np.nan
    data["within_syndrome_significant_fdr"] = False

    for _, indices in data.groupby("syndrome_name", observed=True).groups.items():
        indices = list(indices)
        valid = data.loc[indices, "interaction_p_value"].notna()
        valid_indices = data.loc[indices].index[valid].tolist()

        if not valid_indices:
            continue

        reject, q_values, _, _ = multipletests(
            data.loc[valid_indices, "interaction_p_value"], alpha=FDR_ALPHA, method="fdr_bh"
        )

        data.loc[valid_indices, "within_syndrome_fdr_q_value"] = q_values
        data.loc[valid_indices, "within_syndrome_significant_fdr"] = reject

    return data


def add_global_fdr(data):
    data = data.copy()
    data["global_fdr_q_value"] = np.nan
    data["global_significant_fdr"] = False

    valid_indices = data.index[data["interaction_p_value"].notna()].tolist()

    if valid_indices:
        reject, q_values, _, _ = multipletests(
            data.loc[valid_indices, "interaction_p_value"], alpha=FDR_ALPHA, method="fdr_bh"
        )

        data.loc[valid_indices, "global_fdr_q_value"] = q_values
        data.loc[valid_indices, "global_significant_fdr"] = reject

    return data


# ============================================================
# LOAD DATA
# ============================================================

print("Loading data...")

healthy = pd.read_csv(FAIRFACE_HEALTHY_CSV)
disease = pd.read_csv(GMDB_DISEASE_CSV)

if REQUIRE_FRONTAL:
    if "frontal_ok" in healthy.columns:
        healthy = healthy[healthy["frontal_ok"] == True].copy()

    if "frontal_ok" in disease.columns:
        disease = disease[disease["frontal_ok"] == True].copy()

disease[SYNDROME_COLUMN] = disease[SYNDROME_COLUMN].astype("string").str.strip()

# Harmonized layers
healthy["layer1_group"] = healthy[FAIRFACE_RACE_COLUMN].apply(map_fairface_layer1)
disease["layer1_group"] = disease[GMDB_ETHNICITY_COLUMN].apply(map_gmdb_layer1)

healthy["layer2_group"] = healthy[FAIRFACE_RACE_COLUMN].apply(map_fairface_layer2)
disease["layer2_group"] = [
    map_gmdb_layer2(ethnicity, subcategory)
    for ethnicity, subcategory in zip(disease[GMDB_ETHNICITY_COLUMN], disease[GMDB_SUBCATEGORY_COLUMN])
]

# Exact FaceKit features
feature_columns = list(healthy.loc[:, FIRST_FACEKIT_FEATURE:LAST_FACEKIT_FEATURE].columns)

if len(feature_columns) != 125:
    raise ValueError(f"Expected 125 FaceKit features, found {len(feature_columns)}.")

for feature in feature_columns:
    healthy[feature] = pd.to_numeric(healthy[feature], errors="coerce")
    disease[feature] = pd.to_numeric(disease[feature], errors="coerce")

print("Healthy rows:", len(healthy))
print("Disease rows:", len(disease))
print("FaceKit features:", len(feature_columns))


# ============================================================
# ANALYSIS C
# ============================================================


def run_analysis_c(layer_name, group_column):
    print("\n" + "=" * 60)
    print(layer_name.upper())
    print("=" * 60)

    layer_output = Path(OUTPUT_FOLDER) / layer_name
    layer_output.mkdir(parents=True, exist_ok=True)

    healthy_layer = healthy[healthy[group_column].notna()].copy()
    disease_layer = disease[disease[group_column].notna() & disease[SYNDROME_COLUMN].notna()].copy()

    healthy_counts = healthy_layer[group_column].value_counts()
    eligible_healthy_groups = healthy_counts[healthy_counts >= MIN_HEALTHY_PER_GROUP].index.tolist()

    print("Eligible healthy groups:", eligible_healthy_groups)

    rows = []
    syndromes = sorted(disease_layer[SYNDROME_COLUMN].dropna().unique())

    for syndrome_index, syndrome in enumerate(syndromes, start=1):
        disease_syndrome = disease_layer[disease_layer[SYNDROME_COLUMN] == syndrome].copy()

        disease_counts = disease_syndrome[group_column].value_counts()

        eligible_disease_groups = [
            group for group in eligible_healthy_groups if disease_counts.get(group, 0) >= MIN_DISEASE_PER_GROUP
        ]

        if len(eligible_disease_groups) < 2:
            continue

        group_pairs = list(itertools.combinations(eligible_disease_groups, 2))

        print(
            f"[{syndrome_index}/{len(syndromes)}] {syndrome}: "
            f"{len(eligible_disease_groups)} groups, {len(group_pairs)} pair(s)"
        )

        for group_1, group_2 in group_pairs:
            healthy_pair = healthy_layer[healthy_layer[group_column].isin([group_1, group_2])].copy()

            disease_pair = disease_syndrome[disease_syndrome[group_column].isin([group_1, group_2])].copy()

            for feature in feature_columns:
                h1 = healthy_pair.loc[healthy_pair[group_column] == group_1, feature].dropna().to_numpy()

                h2 = healthy_pair.loc[healthy_pair[group_column] == group_2, feature].dropna().to_numpy()

                d1 = disease_pair.loc[disease_pair[group_column] == group_1, feature].dropna().to_numpy()

                d2 = disease_pair.loc[disease_pair[group_column] == group_2, feature].dropna().to_numpy()

                if (
                    len(h1) < MIN_HEALTHY_PER_GROUP
                    or len(h2) < MIN_HEALTHY_PER_GROUP
                    or len(d1) < MIN_DISEASE_PER_GROUP
                    or len(d2) < MIN_DISEASE_PER_GROUP
                ):
                    continue

                # Build the 2 × 2 interaction dataset.
                interaction_df = pd.DataFrame(
                    {
                        "value": np.concatenate([h1, h2, d1, d2]),
                        "disease": np.concatenate(
                            [np.zeros(len(h1)), np.zeros(len(h2)), np.ones(len(d1)), np.ones(len(d2))]
                        ),
                        "group2": np.concatenate(
                            [np.zeros(len(h1)), np.ones(len(h2)), np.zeros(len(d1)), np.ones(len(d2))]
                        ),
                    }
                )

                try:
                    model = smf.ols("value ~ disease + group2 + disease:group2", data=interaction_df).fit(
                        cov_type="HC3"
                    )

                    interaction_beta = float(model.params["disease:group2"])
                    interaction_se = float(model.bse["disease:group2"])
                    interaction_p = float(model.pvalues["disease:group2"])

                    ci = model.conf_int().loc["disease:group2"]
                    interaction_ci_low = float(ci.iloc[0])
                    interaction_ci_high = float(ci.iloc[1])

                except Exception:
                    interaction_beta = np.nan
                    interaction_se = np.nan
                    interaction_p = np.nan
                    interaction_ci_low = np.nan
                    interaction_ci_high = np.nan

                g1 = hedges_g(d1, h1)
                g2 = hedges_g(d2, h2)
                delta_g = g2 - g1 if np.isfinite(g1) and np.isfinite(g2) else np.nan

                mean_h1 = float(np.mean(h1))
                mean_h2 = float(np.mean(h2))
                mean_d1 = float(np.mean(d1))
                mean_d2 = float(np.mean(d2))

                disease_effect_group1 = mean_d1 - mean_h1
                disease_effect_group2 = mean_d2 - mean_h2

                direction_reversal = bool(
                    np.isfinite(g1)
                    and np.isfinite(g2)
                    and np.sign(g1) != 0
                    and np.sign(g2) != 0
                    and np.sign(g1) != np.sign(g2)
                )

                same_direction = bool(
                    np.isfinite(g1) and np.isfinite(g2) and np.sign(g1) == np.sign(g2) and np.sign(g1) != 0
                )

                rows.append(
                    {
                        "syndrome_name": syndrome,
                        "feature": feature,
                        "group_1": group_1,
                        "group_2": group_2,
                        "n_healthy_group_1": len(h1),
                        "n_healthy_group_2": len(h2),
                        "n_disease_group_1": len(d1),
                        "n_disease_group_2": len(d2),
                        "evidence_level": evidence_level(len(d1), len(d2)),
                        "healthy_mean_group_1": mean_h1,
                        "healthy_mean_group_2": mean_h2,
                        "disease_mean_group_1": mean_d1,
                        "disease_mean_group_2": mean_d2,
                        "raw_disease_effect_group_1": disease_effect_group1,
                        "raw_disease_effect_group_2": disease_effect_group2,
                        "hedges_g_group_1": g1,
                        "hedges_g_group_2": g2,
                        "delta_hedges_g_group2_minus_group1": delta_g,
                        "abs_delta_hedges_g": abs(delta_g) if np.isfinite(delta_g) else np.nan,
                        "interaction_beta_group2_minus_group1": interaction_beta,
                        "interaction_standard_error": interaction_se,
                        "interaction_ci_lower": interaction_ci_low,
                        "interaction_ci_upper": interaction_ci_high,
                        "interaction_ci_excludes_zero": bool(
                            np.isfinite(interaction_ci_low)
                            and np.isfinite(interaction_ci_high)
                            and (interaction_ci_low > 0 or interaction_ci_high < 0)
                        ),
                        "interaction_p_value": interaction_p,
                        "same_effect_direction": same_direction,
                        "direction_reversal": direction_reversal,
                    }
                )

    results = pd.DataFrame(rows)

    if results.empty:
        print("No valid Analysis C comparisons for this layer.")
        return

    # FDR corrections
    results = add_fdr_by_syndrome(results)
    results = add_global_fdr(results)

    # High-priority interaction:
    # - statistically significant within syndrome
    # - robust CI excludes zero
    # - difference in standardized disease effects is sizeable
    # - moderate or primary evidence only

    results["high_priority_interaction"] = (
        results["within_syndrome_significant_fdr"]
        & results["interaction_ci_excludes_zero"]
        & (results["abs_delta_hedges_g"] >= MIN_ABS_DELTA_G_FOR_PRIORITY)
        & results["evidence_level"].isin(["primary", "moderate"])
    )

    # Conservative final flag requiring global FDR too.
    results["high_confidence_interaction"] = (
        results["global_significant_fdr"]
        & results["interaction_ci_excludes_zero"]
        & (results["abs_delta_hedges_g"] >= MIN_ABS_DELTA_G_FOR_PRIORITY)
        & results["evidence_level"].isin(["primary", "moderate"])
    )

    results = results.sort_values(
        ["high_confidence_interaction", "high_priority_interaction", "global_fdr_q_value", "abs_delta_hedges_g"],
        ascending=[False, False, True, False],
    )

    # Full results
    results.to_csv(layer_output / "analysis_C_all_interactions.csv", index=False)

    # Significant within syndrome
    results[results["within_syndrome_significant_fdr"]].to_csv(
        layer_output / "analysis_C_significant_within_syndrome.csv", index=False
    )

    # Global FDR significant
    results[results["global_significant_fdr"]].to_csv(
        layer_output / "analysis_C_global_FDR_significant.csv", index=False
    )

    # Best candidates
    results[results["high_priority_interaction"]].to_csv(
        layer_output / "analysis_C_high_priority_interactions.csv", index=False
    )

    results[results["high_confidence_interaction"]].to_csv(
        layer_output / "analysis_C_high_confidence_interactions.csv", index=False
    )

    # Direction reversals
    results[results["direction_reversal"]].to_csv(layer_output / "analysis_C_direction_reversals.csv", index=False)

    # Syndrome summary
    summary = (
        results.groupby("syndrome_name", observed=True)
        .agg(
            comparisons_tested=("feature", "count"),
            significant_interactions=("within_syndrome_significant_fdr", "sum"),
            global_significant_interactions=("global_significant_fdr", "sum"),
            high_priority_interactions=("high_priority_interaction", "sum"),
            high_confidence_interactions=("high_confidence_interaction", "sum"),
            direction_reversals=("direction_reversal", "sum"),
            max_abs_delta_g=("abs_delta_hedges_g", "max"),
            median_abs_delta_g=("abs_delta_hedges_g", "median"),
        )
        .reset_index()
        .sort_values(
            ["high_confidence_interactions", "high_priority_interactions", "max_abs_delta_g"],
            ascending=[False, False, False],
        )
    )

    summary.to_csv(layer_output / "analysis_C_syndrome_summary.csv", index=False)

    print("\nResults:")
    print("Total interaction tests:", len(results))
    print("Within-syndrome FDR significant:", int(results["within_syndrome_significant_fdr"].sum()))
    print("Global FDR significant:", int(results["global_significant_fdr"].sum()))
    print("High-priority interactions:", int(results["high_priority_interaction"].sum()))
    print("High-confidence interactions:", int(results["high_confidence_interaction"].sum()))
    print("Direction reversals:", int(results["direction_reversal"].sum()))


# ============================================================
# RUN BOTH LAYERS
# ============================================================

Path(OUTPUT_FOLDER).mkdir(parents=True, exist_ok=True)

run_analysis_c(layer_name="layer1_broad_ethnicity", group_column="layer1_group")

run_analysis_c(layer_name="layer2_subcategory", group_column="layer2_group")

print("\nDONE")
print("Saved to:", Path(OUTPUT_FOLDER).resolve())
