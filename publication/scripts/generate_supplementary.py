#!/usr/bin/env python3
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "revision_results/final_20260914"
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT = ROOT / "publication/supplementary"
MANUSCRIPT = ROOT / "manuscript/supplementary"


def clean_syndrome(values): return values.astype(str).str.replace(r"__[0-9a-f]{8}$", "", regex=True)


def truth(values): return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def layer_path(selection, layer, filename):
    return RESULTS / selection / "ethnicity_analysis" / layer / filename


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(EVIDENCE / "feature_region_family_mapping.csv")
    map_cols = ["feature", "anatomical_region", "feature_family", "measurement_description", "normalization_or_formula_context"]
    if len(mapping) != 125 or mapping.feature.nunique() != 125: raise ValueError("Feature mapping must contain the exact 125-feature schema")

    main_results = pd.read_csv(RESULTS / "age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    if len(main_results) != 6125 or main_results.feature.nunique() != 125: raise ValueError("Expected 6,125 tests using 125 features")
    s1 = main_results.merge(mapping[map_cols], on="feature", how="left", validate="many_to_one")
    s1["syndrome"] = clean_syndrome(s1.disease)
    s1["effect_direction"] = np.select([s1.hedges_g.gt(0), s1.hedges_g.lt(0)], ["positive", "negative"], default="zero")
    s1["global_significance_flag"] = s1.mannwhitney_q_fdr_global.lt(.05)
    s1["abs_g_ge_0_8_flag"] = s1.hedges_g.abs().ge(.8)
    s1 = s1.rename(columns={"hedges_g_ci_95_low": "ci_95_lower", "hedges_g_ci_95_high": "ci_95_upper", "mannwhitney_p": "raw_p",
                            "mannwhitney_q_fdr_global": "global_bh_q", "normalization_or_formula_context": "normalization_reference"})
    s1[["syndrome", "feature", "n_disease", "n_healthy", "effect_direction", "hedges_g", "ci_95_lower", "ci_95_upper", "raw_p",
        "global_bh_q", "global_significance_flag", "abs_g_ge_0_8_flag", "anatomical_region", "feature_family", "normalization_reference"]].rename(
        columns={"n_healthy": "n_reference"}).to_csv(OUTPUT / "Supplementary_Table_S1_main_associations.csv", index=False)

    omnibus_rows = []
    for layer, layer_name in [("layer1_broad_ethnicity", "broad"), ("layer2_subcategory", "detailed")]:
        omnibus = pd.read_csv(layer_path("age_3_9", layer, "analysis_A_healthy_kruskal_all_groups.csv.gz"))
        pairwise = pd.read_csv(layer_path("age_3_9", layer, "analysis_A_healthy_pairwise_groups.csv.gz"))
        pairwise = pairwise[truth(pairwise.within_pair_significant_fdr)]
        summaries = {}
        for feature, rows in pairwise.groupby("feature"):
            summaries[feature] = "; ".join(f"{r.group_1} vs {r.group_2}: g={r.hedges_g:.4g}, q={r.within_pair_fdr_q_value:.4g}" for r in rows.itertuples())
        omnibus["demographic_layer"] = layer_name
        omnibus["follow_up_pairwise_results"] = omnibus.feature.map(summaries).fillna("")
        omnibus_rows.append(omnibus)
    s2 = pd.concat(omnibus_rows, ignore_index=True).rename(columns={"kruskal_H": "kruskal_wallis_statistic", "p_value": "raw_p", "fdr_q_value": "corrected_q", "significant_fdr": "significance"})
    s2[["feature", "demographic_layer", "kruskal_wallis_statistic", "raw_p", "corrected_q", "significance", "follow_up_pairwise_results"]].to_csv(
        OUTPUT / "Supplementary_Table_S2_reference_variation.csv", index=False)

    matched, interactions = [], []
    for layer, layer_name in [("layer1_broad_ethnicity", "broad"), ("layer2_subcategory", "detailed")]:
        b = pd.read_csv(layer_path("age_3_9", layer, "analysis_B_matched_disease_vs_healthy.csv.gz"))
        b["demographic_layer"] = layer_name
        matched.append(b)
        c = pd.read_csv(layer_path("age_3_9", layer, "analysis_C_all_interactions.csv.gz"))
        c["demographic_layer"] = layer_name
        interactions.append(c)
    s3 = pd.concat(matched, ignore_index=True).rename(columns={"syndrome_name": "syndrome", "group": "demographic_group", "n_healthy_matched": "n_reference",
                                                       "hedges_g_ci_lower": "ci_95_lower", "hedges_g_ci_upper": "ci_95_upper", "global_fdr_q_value": "global_q",
                                                       "evidence_level": "evidence_tier"})
    s3[["syndrome", "feature", "demographic_layer", "demographic_group", "n_disease", "n_reference", "hedges_g", "ci_95_lower", "ci_95_upper", "global_q", "evidence_tier"]].to_csv(
        OUTPUT / "Supplementary_Table_S3_matched_comparisons.csv", index=False)
    s4 = pd.concat(interactions, ignore_index=True).rename(columns={"syndrome_name": "syndrome", "interaction_standard_error": "hc3_se", "interaction_ci_lower": "ci_95_lower",
                                                            "interaction_ci_upper": "ci_95_upper", "global_fdr_q_value": "global_q", "evidence_level": "evidence_tier"})
    s4[["syndrome", "feature", "demographic_layer", "group_1", "group_2", "n_disease_group_1", "n_disease_group_2", "interaction_beta_group2_minus_group1",
        "hc3_se", "ci_95_lower", "ci_95_upper", "global_q", "hedges_g_group_1", "hedges_g_group_2", "delta_hedges_g_group2_minus_group1", "evidence_tier"]].to_csv(
        OUTPUT / "Supplementary_Table_S4_interactions.csv", index=False)

    alternative = pd.read_csv(RESULTS / "age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    keep = ["disease", "feature", "hedges_g", "mannwhitney_q_fdr_global"]
    s5 = main_results[keep].merge(alternative[keep], on=["disease", "feature"], suffixes=("_primary", "_alternative"), validate="one_to_one")
    s5["syndrome"] = clean_syndrome(s5.disease)
    s5["primary_significance"] = s5.mannwhitney_q_fdr_global_primary.lt(.05)
    s5["alternative_significance"] = s5.mannwhitney_q_fdr_global_alternative.lt(.05)
    s5["direction_agreement"] = np.sign(s5.hedges_g_primary).eq(np.sign(s5.hedges_g_alternative))
    s5["absolute_change_in_g"] = (s5.hedges_g_alternative - s5.hedges_g_primary).abs()
    s5 = s5.rename(columns={"mannwhitney_q_fdr_global_primary": "primary_q", "mannwhitney_q_fdr_global_alternative": "alternative_q",
                            "hedges_g_primary": "primary_g", "hedges_g_alternative": "alternative_g"})
    s5[["syndrome", "feature", "primary_g", "alternative_g", "primary_q", "alternative_q", "primary_significance", "alternative_significance", "direction_agreement", "absolute_change_in_g"]].to_csv(
        OUTPUT / "Supplementary_Table_S5_photo_sensitivity.csv", index=False)

    def status(feature):
        if re.search(r"_asym$", feature): return "bilateral asymmetry"
        if re.search(r"_mean$", feature): return "bilateral aggregate"
        if re.search(r"_l$", feature): return "left-sided"
        if re.search(r"_r$", feature): return "right-sided"
        return "single or derived aggregate"

    s6 = mapping.copy()
    s6["bilateral_aggregate_status"] = s6.feature.map(status)
    s6.rename(columns={"measurement_description": "geometric_definition", "normalization_or_formula_context": "normalization_denominator_reference"})[[
        "feature", "anatomical_region", "feature_family", "bilateral_aggregate_status", "geometric_definition", "normalization_denominator_reference"]].to_csv(
        OUTPUT / "Supplementary_Table_S6_feature_definitions.csv", index=False)
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    for path in OUTPUT.glob("Supplementary_Table_S*.csv"): shutil.copy2(path, MANUSCRIPT / path.name)
    print("Generated S1-S6 supplementary tables")


if __name__ == "__main__": main()
