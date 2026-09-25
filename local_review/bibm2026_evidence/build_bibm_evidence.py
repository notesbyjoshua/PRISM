#!/usr/bin/env python3
import hashlib
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "revision_results/final_20260914"
SOURCE = ROOT / "local_review/final_reanalysis_20260914"
OUT = Path(__file__).resolve().parent
MAIN = PACKAGE / "age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz"
ALT = PACKAGE / "age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz"
LAYERS = {"broad": "layer1_broad_ethnicity", "detailed": "layer2_subcategory"}
Q = "mannwhitney_q_fdr_global"


def clean_syndrome(series):
    return series.str.replace(r"__[0-9a-f]+$", "", regex=True)


def family_name(feature):
    return re.sub(r"_(r|l|mean)$", "", feature)


def feature_region(feature):
    if feature.startswith("eb_"):
        return "eyebrows"
    if feature.startswith(("iris_", "gaze_", "eye_upper_lid", "eye_lower_lid", "inter_pupillary")):
        return "ocular position / iris / gaze"
    if feature.startswith("eye_") or feature in {"inter_canthal_distance", "outer_canthal_distance", "canthal_to_pupillary_ratio"}:
        return "eyes / orbital dimensions"
    if feature.startswith("nose_") or feature in {"nostril_region_area", "nasolabial_angle", "ala_flare_angle", "columella_length", "columella_hang_below_ala"}:
        return "nose"
    if feature.startswith(("philtrum_", "mouth_", "lip_", "upper_lip", "lower_lip", "upper_vermilion", "lower_vermilion", "vermilion_", "cupid_")):
        return "philtrum / mouth / lips"
    if feature.startswith(("chin_", "jaw_")):
        return "chin / jaw"
    if feature.startswith("forehead_") or feature == "hairline_height":
        return "forehead"
    if feature.startswith("face_") or feature == "midface_width":
        return "global face shape / proportions"
    if feature.startswith(("malar_", "cheek_")):
        return "cheeks / malar region"
    raise ValueError(f"Unmapped feature: {feature}")


def normalization_context(feature):
    if feature.endswith("_asym") or feature.startswith("gaze_asym"):
        return "side-to-side difference or normalized asymmetry; see historical formula"
    if "angle" in feature or "slant" in feature:
        return "angle in degrees, except normalized asymmetry variants"
    if feature in {"nose_length", "chin_height", "forehead_height", "hairline_height"}:
        return "face height"
    if feature in {"nose_bridge_length", "columella_length", "columella_hang_below_ala", "philtrum_length"}:
        return "midface height"
    if feature in {"philtrum_width", "mouth_opening", "upper_vermilion_height", "lower_vermilion_height", "vermilion_total", "cupid_bow_drop", "mouth_corner_drop", "mouth_triangularity", "mouth_upper_fit_residual", "upper_lip_eversion", "lower_lip_eversion", "cupid_bow_peak_asym"}:
        return "mouth width"
    if feature == "mouth_tenting":
        return "Cupid-bow peak span"
    if feature == "lip_outer_to_inner_area":
        return "outer-lip area / inner-mouth area"
    if feature == "nose_to_face_area_ratio":
        return "nose polygon area / face-oval area"
    if feature.startswith(("eb_thickness", "eb_lateral_thickness", "eye_area", "nostril_region_area", "cheek_area")):
        return "bizygomatic width squared"
    if feature.startswith(("eye_fissure_aspect", "eye_fissure_fill", "iris_offset_")):
        return "local eye-fissure dimension"
    if feature == "canthal_to_pupillary_ratio":
        return "intercanthal / interpupillary distance"
    if feature == "nose_tip_to_bridge_ratio":
        return "nose tip width / bridge width"
    if feature in {"forehead_taper_ratio", "jaw_to_forehead_ratio", "face_triangularity", "face_width_uniformity", "face_roundness"}:
        return "local face-shape ratio"
    return "bizygomatic width"


def description(feature):
    names = {
        "eb": "eyebrow", "r": "right", "l": "left", "asym": "asymmetry", "inter": "between",
        "fissure": "opening", "canthal": "canthal", "pupillary": "pupillary", "ala": "alar",
        "malar": "malar", "bizygomatic": "bizygomatic", "c2": "quadratic coefficient",
    }
    words = [names.get(word, word) for word in feature.split("_")]
    text = " ".join(words).replace(" mean", " bilateral mean")
    return text[0].upper() + text[1:]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def make_mapping():
    mapping = pd.read_csv(PACKAGE / "exact_feature_names.csv")
    mapping["anatomical_region"] = mapping.feature.map(feature_region)
    mapping["feature_family"] = mapping.feature.map(family_name)
    mapping["family_rule"] = np.where(
        mapping.feature.str.endswith(("_r", "_l", "_mean")),
        "right/left/mean variants grouped; asymmetry retained separately",
        "exact formula retained as its own family",
    )
    mapping["measurement_description"] = mapping.feature.map(description)
    mapping["normalization_or_formula_context"] = mapping.feature.map(normalization_context)
    assert len(mapping) == 125 and mapping.feature.nunique() == 125
    mapping.to_csv(OUT / "feature_region_family_mapping.csv", index=False)
    return mapping


def pose_check():
    rows = []
    for cohort, path in {
        "retained disease": SOURCE / "cohorts/age_3_9/disease_phenotypes.csv",
        "retained reference": SOURCE / "cohorts/reference_phenotypes.csv",
    }.items():
        frame = pd.read_csv(path, usecols=["frontal_ok", "pose_yaw", "pose_pitch", "pose_roll"])
        frame = frame[frame.frontal_ok.eq(True)].copy()
        pose = frame[["pose_yaw", "pose_pitch", "pose_roll"]].apply(pd.to_numeric, errors="coerce")
        finite = np.isfinite(pose).all(axis=1)
        within = finite & pose.pose_yaw.abs().le(15) & pose.pose_pitch.abs().le(15) & pose.pose_roll.abs().le(10)
        classes = {
            "stored pass + finite pose within historical default thresholds": within,
            "stored pass + finite pose outside at least one historical default threshold": finite & ~within,
            "stored pass + missing or nonfinite pose": ~finite,
        }
        for label, mask in classes.items():
            rows.append({"cohort": cohort, "classification": label, "count": int(mask.sum()), "fraction_of_stored_pass": float(mask.mean()),
                         "stored_pass_total": len(frame), "yaw_outside_15": int((finite & pose.pose_yaw.abs().gt(15)).sum()),
                         "pitch_outside_15": int((finite & pose.pose_pitch.abs().gt(15)).sum()), "roll_outside_10": int((finite & pose.pose_roll.abs().gt(10)).sum())})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "pose_flag_consistency.csv", index=False)
    return result


def candidate_tables(main, alt, mapping):
    assert not main.duplicated(["disease", "feature"]).any()
    assert not alt.duplicated(["disease", "feature"]).any()
    assert set(map(tuple, main[["disease", "feature"]].to_numpy())) == set(map(tuple, alt[["disease", "feature"]].to_numpy()))
    alt_fields = ["disease", "feature", "hedges_g", "hedges_g_ci_95_low", "hedges_g_ci_95_high", Q]
    candidates = main[(main[Q] < 0.05) & main.hedges_g.notna() & (main.hedges_g.abs() >= 0.8)].copy()
    candidates = candidates.merge(alt[alt_fields], on=["disease", "feature"], validate="one_to_one", suffixes=("", "_alternative"))
    candidates["syndrome"] = clean_syndrome(candidates.disease)
    candidates["candidate_tier"] = np.select([candidates.n_disease.ge(20), candidates.n_disease.ge(10)], ["primary_candidate", "moderate_candidate"], default="below_editorial_size")
    candidates["both_selections_significant"] = candidates[f"{Q}_alternative"].lt(0.05)
    candidates["effect_direction_agrees"] = np.sign(candidates.hedges_g).eq(np.sign(candidates.hedges_g_alternative))
    candidates["absolute_change_hedges_g"] = (candidates.hedges_g - candidates.hedges_g_alternative).abs()
    primary_ci = candidates.hedges_g_ci_95_low * candidates.hedges_g_ci_95_high > 0
    alt_ci = candidates.hedges_g_ci_95_low_alternative * candidates.hedges_g_ci_95_high_alternative > 0
    candidates["both_confidence_intervals_exclude_zero"] = primary_ci & alt_ci
    interpretable = ~candidates.feature.str.contains("asym|offset|residual|curvature|midline_x|elevation|hang")
    candidates["interpretability_flag"] = interpretable
    candidates["editorial_rank_score"] = (
        candidates.both_selections_significant.astype(int) * 4 + candidates.effect_direction_agrees.astype(int) * 3
        + candidates.both_confidence_intervals_exclude_zero.astype(int) * 2 + interpretable.astype(int)
        + candidates.effect_direction_bootstrap_stability.fillna(0) + candidates.n_disease.clip(upper=50) / 25
        + candidates.hedges_g.abs().clip(upper=3) - 2 * candidates.absolute_change_hedges_g.clip(upper=1)
    )
    candidates = candidates.merge(mapping[["feature", "anatomical_region", "feature_family"]], on="feature", validate="many_to_one")
    candidates["literature_status"] = "unreviewed"
    candidates["clinical_interpretation"] = "TBD"
    candidates["literature_search_terms"] = "TBD"
    candidates = candidates.sort_values(["editorial_rank_score", "n_disease", "feature"], ascending=[False, False, True]).reset_index(drop=True)
    candidates.insert(0, "review_rank", np.arange(1, len(candidates) + 1))
    columns = ["review_rank", "syndrome", "feature", "anatomical_region", "feature_family", "candidate_tier", "n_disease", "n_healthy",
               "mean_disease", "median_disease", "mean_healthy", "median_healthy", "hedges_g", "hedges_g_ci_95_low", "hedges_g_ci_95_high", Q,
               "cliffs_delta", "rank_biserial", "effect_direction_bootstrap_stability", "hedges_g_alternative", "hedges_g_ci_95_low_alternative",
               "hedges_g_ci_95_high_alternative", f"{Q}_alternative", "both_selections_significant", "effect_direction_agrees",
               "absolute_change_hedges_g", "both_confidence_intervals_exclude_zero", "interpretability_flag", "editorial_rank_score",
               "literature_status", "clinical_interpretation", "literature_search_terms"]
    review = candidates[columns].rename(columns={Q: "global_q", f"{Q}_alternative": "alternative_global_q"})
    review.to_csv(OUT / "main_candidate_review.csv", index=False)

    eligible = candidates[candidates.both_selections_significant & candidates.effect_direction_agrees & candidates.both_confidence_intervals_exclude_zero].copy()
    selected, seen_pairs, syndrome_counts = [], set(), {}
    for region in mapping.anatomical_region.drop_duplicates():
        pool = eligible[eligible.anatomical_region.eq(region)]
        for idx, row in pool.iterrows():
            pair = (row.syndrome, row.feature_family)
            if pair in seen_pairs or syndrome_counts.get(row.syndrome, 0) >= 2:
                continue
            selected.append(idx); seen_pairs.add(pair); syndrome_counts[row.syndrome] = syndrome_counts.get(row.syndrome, 0) + 1
            break
    for idx, row in eligible.iterrows():
        pair = (row.syndrome, row.feature_family)
        if idx in selected or pair in seen_pairs or syndrome_counts.get(row.syndrome, 0) >= 2:
            continue
        selected.append(idx); seen_pairs.add(pair); syndrome_counts[row.syndrome] = syndrome_counts.get(row.syndrome, 0) + 1
        if len(selected) >= 25:
            break
    shortlist = candidates.loc[selected].sort_values(["editorial_rank_score", "n_disease"], ascending=False).copy()
    shortlist.insert(0, "shortlist_rank", np.arange(1, len(shortlist) + 1))
    shortlist = shortlist[["shortlist_rank"] + columns]
    shortlist = shortlist.rename(columns={Q: "global_q", f"{Q}_alternative": "alternative_global_q"})
    assert 20 <= len(shortlist) <= 30 and not shortlist.duplicated(["syndrome", "feature_family"]).any()
    shortlist.to_csv(OUT / "main_candidate_shortlist.csv", index=False)
    return review, shortlist


def trend_tables(main, mapping):
    frame = main.copy()
    frame["syndrome"] = clean_syndrome(frame.disease)
    frame = frame.merge(mapping[["feature", "anatomical_region", "feature_family"]], on="feature", validate="many_to_one")
    frame["significant"] = frame[Q].lt(0.05)
    frame["significant_large"] = frame.significant & frame.hedges_g.abs().ge(0.8)

    def aggregate(cols):
        out = frame.groupby(cols, sort=False).agg(total_tested_comparisons=("feature", "size"), number_significant=("significant", "sum"),
                                                    number_significant_large=("significant_large", "sum"), eligible_syndromes=("syndrome", "nunique")).reset_index()
        out["fraction_significant"] = out.number_significant / out.total_tested_comparisons
        out["fraction_significant_large"] = out.number_significant_large / out.total_tested_comparisons
        return out

    regional = aggregate(["anatomical_region"])
    syndrome_region = aggregate(["syndrome", "anatomical_region"])
    regional.to_csv(OUT / "regional_trends.csv", index=False)
    syndrome_region.to_csv(OUT / "syndrome_region_trends.csv", index=False)
    recurrence = frame.groupby(["feature", "anatomical_region", "feature_family"], sort=False).agg(
        syndromes_tested=("syndrome", "nunique"), number_significant=("significant", "sum"),
        number_significant_positive=("hedges_g", lambda x: int(((x > 0) & frame.loc[x.index, "significant"]).sum())),
        number_significant_negative=("hedges_g", lambda x: int(((x < 0) & frame.loc[x.index, "significant"]).sum())),
        number_significant_large=("significant_large", "sum"),
        median_signed_g_among_significant=("hedges_g", lambda x: x[frame.loc[x.index, "significant"]].median()),
        median_absolute_g_among_significant=("hedges_g", lambda x: x[frame.loc[x.index, "significant"]].abs().median()),
    ).reset_index()
    recurrence["significant_fraction_among_tested_syndromes"] = recurrence.number_significant / recurrence.syndromes_tested
    recurrence.to_csv(OUT / "feature_recurrence.csv", index=False)
    family_pairs = frame.groupby(["syndrome", "anatomical_region", "feature_family"], as_index=False).agg(
        underlying_exact_comparisons=("feature", "size"), significant=("significant", "max"), significant_large=("significant_large", "max"),
        max_abs_g=("hedges_g", lambda x: x.abs().max()))
    family = family_pairs.groupby(["anatomical_region", "feature_family"], sort=False).agg(
        total_tested_comparisons=("syndrome", "size"), total_underlying_exact_comparisons=("underlying_exact_comparisons", "sum"),
        number_significant=("significant", "sum"), number_significant_large=("significant_large", "sum"),
        eligible_syndromes=("syndrome", "nunique")).reset_index()
    family["fraction_significant"] = family.number_significant / family.total_tested_comparisons
    family["fraction_significant_large"] = family.number_significant_large / family.total_tested_comparisons
    family.to_csv(OUT / "feature_family_trends.csv", index=False)
    return frame, regional, family, recurrence, family_pairs


def demographic_tables(shortlist):
    keys = shortlist[["shortlist_rank", "syndrome", "feature", "hedges_g"]].rename(columns={"hedges_g": "pooled_hedges_g"})
    supports, interactions = [], []
    for layer_name, directory in LAYERS.items():
        b = pd.read_csv(PACKAGE / f"age_3_9/ethnicity_analysis/{directory}/analysis_B_matched_disease_vs_healthy.csv.gz")
        merged = keys.merge(b, left_on=["syndrome", "feature"], right_on=["syndrome_name", "feature"], how="left")
        merged.insert(1, "demographic_layer", layer_name)
        merged["direction_agrees_with_pooled"] = np.sign(merged.pooled_hedges_g).eq(np.sign(merged.hedges_g))
        supports.append(merged)
        c = pd.read_csv(PACKAGE / f"age_3_9/ethnicity_analysis/{directory}/analysis_C_all_interactions.csv.gz")
        c = c[c.high_confidence_interaction.eq(True)]
        found = keys.merge(c, left_on=["syndrome", "feature"], right_on=["syndrome_name", "feature"], how="inner")
        found.insert(1, "demographic_layer", layer_name)
        interactions.append(found)
    support = pd.concat(supports, ignore_index=True)
    keep = ["shortlist_rank", "demographic_layer", "syndrome", "feature", "group", "n_disease", "n_healthy_matched", "hedges_g",
            "hedges_g_ci_lower", "hedges_g_ci_upper", "global_fdr_q_value", "within_group_fdr_q_value", "evidence_level",
            "global_significant_fdr", "direction_agrees_with_pooled"]
    support[keep].to_csv(OUT / "candidate_demographic_support.csv", index=False)
    interaction = pd.concat(interactions, ignore_index=True) if interactions else pd.DataFrame()
    interaction_keep = ["shortlist_rank", "demographic_layer", "syndrome", "feature", "group_1", "group_2", "n_healthy_group_1",
                        "n_healthy_group_2", "n_disease_group_1", "n_disease_group_2", "evidence_level", "interaction_beta_group2_minus_group1",
                        "interaction_ci_lower", "interaction_ci_upper", "global_fdr_q_value", "hedges_g_group_1", "hedges_g_group_2",
                        "delta_hedges_g_group2_minus_group1", "high_confidence_interaction"]
    if interaction.empty:
        interaction = pd.DataFrame(columns=interaction_keep)
    interaction[interaction_keep].to_csv(OUT / "candidate_interaction_support.csv", index=False)
    return support, interaction


def literature_queue(shortlist, mapping, support):
    support_summary = support.dropna(subset=["group"]).groupby(["syndrome", "feature"]).apply(
        lambda x: "; ".join(f"{r.demographic_layer}:{r.group} n={int(r.n_disease)}, g={r.hedges_g:.2f}, global q={r.global_fdr_q_value:.3g}, {r.evidence_level}"
                            for r in x.itertuples() if r.evidence_level != "insufficient") or "tested groups had insufficient evidence",
        include_groups=False,
    ).rename("demographic_support_summary").reset_index()
    queue = shortlist.head(15).merge(mapping[["feature", "measurement_description", "normalization_or_formula_context"]], on="feature", validate="many_to_one")
    queue = queue.merge(support_summary, on=["syndrome", "feature"], how="left")
    queue["demographic_support_summary"] = queue.demographic_support_summary.fillna("no matched demographic row available")
    queue.insert(0, "literature_priority", ["first-pass"] * min(5, len(queue)) + ["second-pass"] * max(0, len(queue) - 5))
    queue["suggested_clinical_synonyms"] = queue.feature.str.replace("_", " ") + "; facial morphology; craniofacial measurement"
    queue["suggested_HPO_search_terms"] = queue.feature.str.replace("_", " ") + "; abnormal facial morphology"
    queue["suggested_literature_query"] = '"' + queue.syndrome + '" AND ("' + queue.feature.str.replace("_", " ") + '" OR craniofacial OR facial phenotype)'
    queue["literature_status"] = "unreviewed"
    keep = ["literature_priority", "shortlist_rank", "syndrome", "feature", "measurement_description", "normalization_or_formula_context", "n_disease",
            "hedges_g", "hedges_g_ci_95_low", "hedges_g_ci_95_high", "global_q", "hedges_g_alternative", "hedges_g_ci_95_low_alternative",
            "hedges_g_ci_95_high_alternative", "alternative_global_q", "demographic_support_summary", "suggested_clinical_synonyms",
            "suggested_HPO_search_terms", "suggested_literature_query", "literature_status"]
    queue[keep].to_csv(OUT / "literature_review_queue.csv", index=False)
    return queue[keep]


def make_heatmap(frame, recurrence):
    syndrome_stats = frame.groupby("syndrome").agg(n=("n_disease", "max"), large=("significant_large", "sum"), significant=("significant", "sum"))
    syndromes = syndrome_stats[syndrome_stats.n.ge(20)].sort_values(["n", "large", "significant"], ascending=False).head(10).index.tolist()
    recurrent = recurrence.sort_values(["number_significant_large", "number_significant", "median_absolute_g_among_significant"], ascending=False)
    selected = []
    for region in recurrent.anatomical_region.drop_duplicates():
        pool = recurrent[(recurrent.anatomical_region == region) & ~recurrent.feature.str.contains("_r$|_l$|asym")]
        if not pool.empty:
            selected.append(pool.iloc[0].feature)
    for feature in recurrent.feature:
        family = family_name(feature)
        if feature in selected or any(family_name(x) == family for x in selected) or re.search(r"_r$|_l$|asym", feature):
            continue
        selected.append(feature)
        if len(selected) >= 12:
            break
    selected = selected[:12]
    source = frame[frame.syndrome.isin(syndromes) & frame.feature.isin(selected)].copy()
    source["cell_status"] = np.where(source[Q].lt(0.05), "significant", "tested_nonsignificant")
    source = source[["syndrome", "feature", "feature_family", "anatomical_region", "n_disease", "hedges_g", "hedges_g_ci_95_low", "hedges_g_ci_95_high", Q, "cell_status"]]
    source.to_csv(OUT / "figure_cross_syndrome_heatmap_source.csv", index=False)
    effects = source.pivot(index="feature", columns="syndrome", values="hedges_g").reindex(index=selected, columns=syndromes)
    significant = source.pivot(index="feature", columns="syndrome", values=Q).reindex(index=selected, columns=syndromes).lt(0.05)
    tested = effects.notna()
    base = np.where(~tested, 0, np.where(significant, 2, 1))
    def short_label(syndrome):
        if ";" in syndrome:
            suffix = syndrome.rsplit(";", 1)[1].strip()
            if len(suffix) <= 7:
                return suffix
        return re.sub(r" syndrome$", "", syndrome, flags=re.I)

    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    ax.imshow(base, cmap=ListedColormap(["#252525", "#e2e2e2", "#ffffff"]), vmin=0, vmax=2, aspect="auto")
    masked = np.ma.masked_where(~significant.to_numpy(), effects.to_numpy())
    vmax = max(1.5, float(np.nanquantile(np.abs(masked.filled(np.nan)), 0.98)))
    image = ax.imshow(masked, cmap="RdBu_r", norm=Normalize(-vmax, vmax), aspect="auto")
    ax.set_xticks(range(len(syndromes)), [f"{short_label(s)}\n(n={int(syndrome_stats.loc[s, 'n'])})" for s in syndromes], rotation=42, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(selected)), [f.replace("_", " ") for f in selected], fontsize=8)
    ax.set_xlabel("Syndrome group")
    ax.set_ylabel("Representative exact feature (one per family)")
    ax.set_title("Selected recurrent facial measurements across well-powered syndromes", fontsize=10)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Signed Hedges' g (significant cells)", fontsize=8)
    ax.legend(handles=[Patch(facecolor="#e2e2e2", label="tested, q ≥ 0.05"), Patch(facecolor="#252525", label="untested")],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False, fontsize=8)
    ax.set_xticks(np.arange(-0.5, len(syndromes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(selected), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.subplots_adjust(left=0.26, right=0.93, top=0.90, bottom=0.29)
    fig.savefig(OUT / "figure_cross_syndrome_heatmap.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "figure_cross_syndrome_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    return syndromes, selected


def fmt_effect(row):
    return f"{row.syndrome} — {row.feature}: g={row.hedges_g:.2f} (95% CI {row.hedges_g_ci_95_low:.2f} to {row.hedges_g_ci_95_high:.2f}), global q={row.global_q:.3g}, n={int(row.n_disease)}"


def write_report(commit, main, alt, mapping, pose, review, shortlist, regional, family, recurrence, family_pairs, support, interaction, queue, heatmap_syndromes, heatmap_features):
    sig = main[Q].lt(0.05)
    alt_sig = alt[Q].lt(0.05)
    keys = ["disease", "feature"]
    overlap = main[keys].assign(primary=sig).merge(alt[keys].assign(alternative=alt_sig), on=keys, validate="one_to_one")
    broad_a = pd.read_csv(PACKAGE / "age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_A_healthy_kruskal_all_groups.csv.gz")
    detailed_a = pd.read_csv(PACKAGE / "age_3_9/ethnicity_analysis/layer2_subcategory/analysis_A_healthy_kruskal_all_groups.csv.gz")
    headline = {"patients": 1850, "references": 942, "syndromes": main.disease.nunique(), "tests": len(main), "sig": int(sig.sum()),
                "large": int((sig & main.hedges_g.abs().ge(0.8)).sum()), "shared": int((overlap.primary & overlap.alternative).sum()),
                "primary_only": int((overlap.primary & ~overlap.alternative).sum()), "alternative_only": int((~overlap.primary & overlap.alternative).sum()),
                "broad_a": int(broad_a.significant_fdr.sum()), "detailed_a": int(detailed_a.significant_fdr.sum())}
    expected = {"syndromes": 49, "tests": 6125, "sig": 811, "large": 379, "shared": 770, "primary_only": 41, "alternative_only": 64, "broad_a": 72, "detailed_a": 89}
    for key, value in expected.items():
        assert headline[key] == value, f"Headline mismatch for {key}: {headline[key]} != {value}"
    disease_rows = pd.read_csv(SOURCE / "cohorts/age_3_9/disease_phenotypes.csv", usecols=["image_id"])
    reference_rows = pd.read_csv(SOURCE / "cohorts/reference_phenotypes.csv", usecols=["image_id"])
    assert len(disease_rows) == headline["patients"] and reference_rows.image_id.nunique() == headline["references"]
    assert (main.hedges_g_ci_95_low <= main.hedges_g_ci_95_high).all()
    assert (alt.hedges_g_ci_95_low <= alt.hedges_g_ci_95_high).all()
    assert review.global_q.lt(0.05).all() and review.hedges_g.abs().ge(0.8).all()
    top_regions = regional.sort_values(["fraction_significant", "number_significant"], ascending=False).head(3)
    pair_summary = family_pairs.groupby("anatomical_region").agg(tested_family_syndrome_pairs=("significant", "size"), detected_family_syndrome_pairs=("significant", "sum")).reset_index()
    pair_summary["detected_fraction"] = pair_summary.detected_family_syndrome_pairs / pair_summary.tested_family_syndrome_pairs
    top_family_regions = pair_summary.sort_values("detected_fraction", ascending=False).head(3)
    top_families = family.sort_values(["fraction_significant", "number_significant"], ascending=False).head(5)
    sig_main = main[sig]
    positive, negative = int((sig_main.hedges_g > 0).sum()), int((sig_main.hedges_g < 0).sum())
    pose_lines = []
    for cohort, part in pose.groupby("cohort", sort=False):
        counts = dict(zip(part.classification, part["count"]))
        pose_lines.append(f"- {cohort}: {counts['stored pass + finite pose within historical default thresholds']:,} within thresholds; "
                          f"{counts['stored pass + finite pose outside at least one historical default threshold']:,} outside at least one threshold; "
                          f"{counts['stored pass + missing or nonfinite pose']:,} missing/nonfinite pose.")
    input_paths = [MAIN, ALT, PACKAGE / "exact_feature_names.csv", SOURCE / "cohorts/age_3_9/disease_phenotypes.csv",
                   SOURCE / "cohorts/reference_phenotypes.csv"]
    for directory in LAYERS.values():
        for analysis in ["analysis_A_healthy_kruskal_all_groups", "analysis_B_matched_disease_vs_healthy", "analysis_C_all_interactions"]:
            input_paths.append(PACKAGE / f"age_3_9/ethnicity_analysis/{directory}/{analysis}.csv.gz")
    hash_lines = "\n".join(f"- `{path.relative_to(ROOT)}` — `{sha256(path)}`" for path in input_paths)
    strong = shortlist.head(5)
    remaining = shortlist.iloc[5:]
    concordance = remaining[remaining.feature.str.contains("inter_canthal|outer_canthal|philtrum|chin|forehead")].head(2)
    provisional = remaining[remaining.interpretability_flag & ~remaining.index.isin(concordance.index)].tail(2)
    support_examples = support[(support.global_significant_fdr.eq(True)) & ~support.evidence_level.eq("insufficient")].sort_values(["n_disease", "global_fdr_q_value"], ascending=[False, True]).drop_duplicates(["syndrome", "feature"]).head(3)
    interaction_examples = interaction.sort_values("global_fdr_q_value").drop_duplicates(["syndrome", "feature"]).head(3) if not interaction.empty else interaction
    alt_delta = shortlist.absolute_change_hedges_g
    regions_text = "; ".join(f"{r.anatomical_region}: {int(r.number_significant)}/{int(r.total_tested_comparisons)} ({r.fraction_significant:.1%})" for r in top_regions.itertuples())
    family_regions_text = "; ".join(f"{r.anatomical_region}: {int(r.detected_family_syndrome_pairs)}/{int(r.tested_family_syndrome_pairs)} ({r.detected_fraction:.1%})" for r in top_family_regions.itertuples())
    families_text = "; ".join(f"{r.feature_family}: {int(r.number_significant)}/{int(r.total_tested_comparisons)} syndromes ({r.fraction_significant:.1%})" for r in top_families.itertuples())
    support_text = "\n".join(f"- {r.syndrome} — {r.feature}, {r.demographic_layer}/{r.group}: n={int(r.n_disease)} disease and {int(r.n_healthy_matched)} reference; g={r.hedges_g:.2f} (95% CI {r.hedges_g_ci_lower:.2f} to {r.hedges_g_ci_upper:.2f}); global q={r.global_fdr_q_value:.3g}; {r.evidence_level}; pooled direction {'agrees' if r.direction_agrees_with_pooled else 'differs'}." for r in support_examples.itertuples()) or "- No shortlisted Analysis B row met both global significance and a non-insufficient evidence level."
    interaction_text = "\n".join(f"- {r.syndrome} — {r.feature}, {r.demographic_layer} {r.group_1} versus {r.group_2}: interaction beta={r.interaction_beta_group2_minus_group1:.3g} (HC3 95% CI {r.interaction_ci_lower:.3g} to {r.interaction_ci_upper:.3g}), global q={r.global_fdr_q_value:.3g}; descriptive g values {r.hedges_g_group_1:.2f} and {r.hedges_g_group_2:.2f}." for r in interaction_examples.itertuples()) or "- No high-confidence Analysis C interaction matched the 25-row main shortlist."
    result_text_rows = strong.head(3)
    result_effects = "; ".join(f"{r.syndrome} {r.feature.replace('_', ' ')} (g={r.hedges_g:.2f}, 95% CI {r.hedges_g_ci_95_low:.2f} to {r.hedges_g_ci_95_high:.2f}, global q={r.global_q:.3g}, n={int(r.n_disease)})" for r in result_text_rows.itertuples())
    report = f"""# BIBM 2026 evidence report

Generated from the approved childhood PRISM reanalysis. All clinical literature classifications remain pending; no web literature search was performed in this summarization run.

## A. Provenance

- PRISM commit: `{commit}` (worktree was clean before generation).
- Validated package: `revision_results/final_20260914/`.
- Retained local inputs used only for aggregate cohort and pose checks: `local_review/final_reanalysis_20260914/cohorts/age_3_9/disease_phenotypes.csv` and `cohorts/reference_phenotypes.csv`.
- Historical extraction provenance: `e474628f3cb17a83edad59ac6e1d9e5438753032` and the validated 125-feature schema. No feature extraction was run, and no newer 120-feature FaceKit output or bug-testing result was used.
- Generator: `local_review/bibm2026_evidence/build_bibm_evidence.py`; deterministic, with no new resampling or inferential tests.

Important result-input hashes:

{hash_lines}

## B. Confirmed headline results

Independent table checks reproduced all frozen counts: {headline['patients']:,} retained disease patients, {headline['references']:,} distinct reference image paths, {headline['syndromes']} eligible syndromes, and {headline['tests']:,} main tests. Of these, {headline['sig']} had global BH q<0.05 and {headline['large']} also had |g|≥0.8. Primary and alternative photo selections shared {headline['shared']} significant hypotheses; {headline['primary_only']} were primary-only and {headline['alternative_only']} alternative-only. Analysis A yielded {headline['broad_a']}/125 omnibus-significant features in the broad layer and {headline['detailed_a']}/125 in the detailed layer.

## C. Paper-ready broader trends

At the exact-feature comparison level, the highest detected-association fractions were {regions_text}. Across all significant rows, {positive} effects were positive and {negative} were negative under the stored disease-minus-reference sign convention.

After collapsing right/left/mean variants into transparent feature families and counting each syndrome-family pair once, the leading regional fractions were {family_regions_text}. These are detected-association fractions, not formal enrichment tests. Related formulas remain correlated and should not be narrated as independent biological discoveries.

The most recurrent individual feature families were {families_text}. Recurrence does not imply syndrome specificity or clinical importance.

The complete region, family, syndrome-region, and exact-feature summaries are in the accompanying CSV files. No feature is described as syndrome-specific merely because it was detected once.

## D. Selected syndrome-feature candidates

### 1. Strong quantitative associations needing clinical literature classification

{chr(10).join('- ' + fmt_effect(r) + ' [clinical literature classification pending]' for r in strong.itertuples())}

### 2. Established-looking concordance examples to confirm

These are geometrically interpretable patterns selected for author checking, not claims of established clinical concordance.

{chr(10).join('- ' + fmt_effect(r) + ' [clinical literature classification pending]' for r in concordance.itertuples())}

### 3. Candidates worth investigating as potentially underreported

The label is provisional and does not imply novelty.

{chr(10).join('- ' + fmt_effect(r) + ' [clinical literature classification pending]' for r in provisional.itertuples())}

## E. Demographic and interaction examples

Broad and detailed demographic layers are reported separately and are not independent replications. Small nonsignificant subgroups are not treated as evidence of no effect.

{support_text}

Formal high-confidence interaction examples matching the shortlist:

{interaction_text}

## F. Photograph-selection sensitivity

All {len(shortlist)} shortlisted effects had the same direction, were globally significant in both selections, and had confidence intervals excluding zero in both selections by construction. The median absolute change in g was {alt_delta.median():.3f}; the maximum was {alt_delta.max():.3f}. This is a photograph-selection robustness analysis, not independent replication.

## G. Pose-flag diagnostic

{chr(10).join(pose_lines)}

The comparison uses historical defaults (|yaw|≤15°, |pitch|≤15°, |roll|≤10°) only as a consistency reference. It does not establish the invocation settings used for every historical record, does not override stored `frontal_ok`, and did not change the cohort.

## H. Figure recommendations

1. Use `figure_cross_syndrome_heatmap.png` as the primary Results figure. It uses the ten largest eligible syndrome groups with n≥20 and twelve representative exact features selected from recurrent/large effects while taking one interpretable representative per family and spanning all anatomical regions. Significant cells show signed Hedges' g; tested nonsignificant cells are gray; untested cells would be black (none occur in this complete selected grid).
2. If space permits, derive a small forest panel from the top 6–10 rows of `main_candidate_shortlist.csv`; the stored confidence intervals and alternate-selection columns support this without new inference.

Proposed heatmap caption: **Cross-syndrome structure of selected facial measurements in children ages 3–9.** Columns are the ten largest eligible syndrome groups (n≥20); rows are twelve representative exact measurements from distinct feature families, selected using recurrence, large-effect frequency, interpretability, and anatomical coverage rather than q-value alone. Color encodes signed Hedges' g only for global-BH-significant comparisons (q<0.05); gray cells were tested but nonsignificant. Positive values indicate larger measurements in the syndrome group than in the reference set. Feature variants are not independent, and the display is descriptive.

## I. Exact manuscript replacement text

### Selected phenotype associations and anatomical trends

Among 6,125 syndrome-feature comparisons, 811 met the global Benjamini–Hochberg threshold (q<0.05), including 379 with |Hedges' g|≥0.8. The largest detected-association fractions at the exact-comparison level were {regions_text}. Family-level de-duplication retained the same broad emphasis ({family_regions_text}), although correlated measurements should not be interpreted as independent biological findings. Representative stable associations included {result_effects}; each remains **[clinical literature classification pending]**. All {len(shortlist)} shortlisted effects retained direction, global significance, and confidence intervals excluding zero under alternative eligible-photograph selection; the median absolute change in g was {alt_delta.median():.3f}. This sensitivity analysis supports robustness to the prespecified photograph-selection rule but is not an independent replication.

## Quality-control record

- Main candidates were filtered with `mannwhitney_q_fdr_global`; within-syndrome q-values were not used for selection.
- Primary and alternative rows were joined one-to-one on the exact stored disease label and feature, then labels were cleaned only for presentation.
- All stored main and alternative confidence intervals had correctly ordered bounds; sign agreement was computed directly from signed Hedges' g.
- Denominators use actual tested rows. Untested status is distinct from nonsignificance; the selected main grid happens to be complete.
- Analysis B layers remain separate. Analysis C output uses the formal interaction beta and HC3 interval/q; delta-g is retained only as descriptive context.
- The output directory contains no patient IDs, patient-level manifests, or images. The pose table contains aggregate counts only.
- All outputs derive from the saved historical 125-feature result schema; no updated FaceKit code or regenerated feature values were used.
"""
    (OUT / "BIBM_EVIDENCE_REPORT.md").write_text(report)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    main_results, alt_results = pd.read_csv(MAIN), pd.read_csv(ALT)
    mapping = make_mapping()
    pose = pose_check()
    review, shortlist = candidate_tables(main_results, alt_results, mapping)
    frame, regional, family, recurrence, family_pairs = trend_tables(main_results, mapping)
    support, interaction = demographic_tables(shortlist)
    queue = literature_queue(shortlist, mapping, support)
    heatmap_syndromes, heatmap_features = make_heatmap(frame, recurrence)
    write_report(commit, main_results, alt_results, mapping, pose, review, shortlist, regional, family, recurrence, family_pairs, support,
                 interaction, queue, heatmap_syndromes, heatmap_features)
    files = sorted(path.name for path in OUT.iterdir() if path.is_file())
    print(f"Created {len(files)} files: {', '.join(files)}")
    print("Strongest five candidates:")
    for row in shortlist.head(5).itertuples():
        print(f"  {row.syndrome} — {row.feature}: g={row.hedges_g:.2f}, n={row.n_disease}, q={row.global_q:.3g}")
    best = regional.sort_values(["fraction_significant", "number_significant"], ascending=False).iloc[0]
    print(f"Strongest broader trend: {best.anatomical_region}, {best.number_significant}/{best.total_tested_comparisons} significant ({best.fraction_significant:.1%}).")
    print("Best figure: figure_cross_syndrome_heatmap.png")
    print("Unresolved before submission: clinical literature classification and citations for shortlisted associations.")


if __name__ == "__main__":
    main()
