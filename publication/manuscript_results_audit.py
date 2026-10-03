#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "revision_results/final_20260914"
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT_CSV = ROOT / "publication/manuscript_results_audit.csv"
OUTPUT_MD = ROOT / "publication/MANUSCRIPT_RESULTS_AUDIT.md"


def truth(values):
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def clean_syndrome(value):
    return re.sub(r"__[0-9a-f]{8}$", "", str(value))


def fmt(value):
    if isinstance(value, (float, np.floating)):
        if math.isnan(value): return "NaN"
        if value and (abs(value) < 1e-4 or abs(value) >= 1e6): return f"{value:.6g}"
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


class Audit:
    def __init__(self): self.rows = []

    def add(self, claim, expected, computed, source, calculation, tolerance=0, notes=""):
        if computed is None or (isinstance(computed, float) and math.isnan(computed)):
            status = "UNRESOLVED"
        elif isinstance(expected, (int, float)) and isinstance(computed, (int, float, np.number)):
            status = "PASS" if abs(float(computed) - float(expected)) <= tolerance else "FAIL"
        else:
            status = "PASS" if str(computed) == str(expected) else "FAIL"
        self.rows.append({"claim_name": claim, "expected_manuscript_value": fmt(expected), "computed_value": fmt(computed),
                          "source_file": source, "calculation": calculation, "status": status, "notes": notes})

    def write(self):
        frame = pd.DataFrame(self.rows)
        frame.to_csv(OUTPUT_CSV, index=False)
        counts = frame.status.value_counts().to_dict()
        lines = ["# Manuscript results audit", "", "This audit reconstructs manuscript claims from the corrected historical 125-feature analysis. It does not use newer FaceKit outputs.", "",
                 f"**Result:** {counts.get('PASS', 0)} PASS, {counts.get('FAIL', 0)} FAIL, {counts.get('UNRESOLVED', 0)} UNRESOLVED.", "",
                 "| Claim | Expected | Computed | Status | Source and calculation |", "|---|---:|---:|:---:|---|"]
        for row in self.rows:
            detail = f"`{row['source_file']}`; {row['calculation']}"
            if row["notes"]: detail += f" Notes: {row['notes']}"
            detail = detail.replace("|", "\\|")
            lines.append(f"| {row['claim_name']} | {row['expected_manuscript_value']} | {row['computed_value']} | {row['status']} | {detail} |")
        lines += ["", "The script exits nonzero if any checked claim fails or cannot be reconstructed.", ""]
        OUTPUT_MD.write_text("\n".join(lines))
        return frame


def row_for(frame, syndrome, feature):
    names = frame["disease"].map(clean_syndrome)
    rows = frame[names.eq(syndrome) & frame.feature.eq(feature)]
    if len(rows) != 1: raise ValueError(f"Expected one result for {syndrome} / {feature}, found {len(rows)}")
    return rows.iloc[0]


def subgroup_row(frame, syndrome, feature, group):
    rows = frame[frame.syndrome_name.eq(syndrome) & frame.feature.eq(feature) & frame.group.eq(group)]
    if len(rows) != 1: raise ValueError(f"Expected one subgroup result for {syndrome} / {feature} / {group}, found {len(rows)}")
    return rows.iloc[0]


def main():
    audit = Audit()
    disease = pd.read_csv(ROOT / "data/master_dataset_disease.csv", low_memory=False)
    reference = pd.read_csv(ROOT / "data/master_dataset_healthy.csv", low_memory=False)
    features = pd.read_csv(RESULTS / "exact_feature_names.csv").feature.tolist()
    main_results = pd.read_csv(RESULTS / "age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    alternative = pd.read_csv(RESULTS / "age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    retention = pd.read_csv(RESULTS / "syndrome_retention_and_main_tests.csv")
    primary_retention = retention[retention.cohort.eq("age_3_9")]
    quality = truth(disease.frontal_ok)
    year, month = pd.to_numeric(disease.age_year, errors="coerce"), pd.to_numeric(disease.age_month, errors="coerce")
    valid = year.notna() & month.notna() & year.ge(0) & month.between(0, 11) & year.mod(1).eq(0) & month.mod(1).eq(0)
    age = (year + month / 12).where(valid & ~(year.eq(0) & month.eq(0)))
    age_eligible = quality & age.ge(3) & age.lt(10)

    audit.add("Disease-image records", 10942, len(disease), "data/master_dataset_disease.csv", "row count")
    audit.add("Stored frontal-quality passes", 9411, int(quality.sum()), "data/master_dataset_disease.csv", "count frontal_ok values parsed as true")
    audit.add("Age-eligible frontal-pass photographs", 2101, int(age_eligible.sum()), "data/master_dataset_disease.csv", "frontal pass and 3 <= age_year + age_month/12 < 10 after validating integer components")
    eligible_patients = disease.loc[age_eligible, "patient_id"].astype("string").str.strip().nunique()
    audit.add("Repeat-patient photographs removed", 251, int(age_eligible.sum() - eligible_patients), "data/master_dataset_disease.csv", "eligible photograph rows minus distinct known patient IDs")
    audit.add("Retained disease patients", 1850, eligible_patients, "data/master_dataset_disease.csv", "distinct patient IDs among age-eligible frontal-pass rows")
    audit.add("Retained syndrome labels", 379, primary_retention.syndrome.nunique(), "revision_results/final_20260914/syndrome_retention_and_main_tests.csv", "distinct syndrome labels in primary cohort rows")
    audit.add("Main-analysis eligible syndrome groups", 49, int(truth(primary_retention.main_eligible).sum()), "revision_results/final_20260914/syndrome_retention_and_main_tests.csv", "count main_eligible=true in primary cohort")
    audit.add("Distinct reference image paths", 942, reference.file.astype("string").str.replace("\\", "/", regex=False).nunique(), "data/master_dataset_healthy.csv", "normalized distinct file paths")
    audit.add("Measurements per eligible syndrome", 125, main_results.feature.nunique(), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "distinct feature names")
    audit.add("Primary syndrome-feature tests", 6125, len(main_results), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "row count")
    audit.add("Historical feature schema", 125, len(features), "revision_results/final_20260914/exact_feature_names.csv", "row count", notes="Exact schema matches all 125 feature columns used by the result table.")

    sig = main_results.mannwhitney_q_fdr_global.lt(.05)
    large = main_results.hedges_g.abs().ge(.8)
    audit.add("Globally significant comparisons", 811, int(sig.sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count global BH q < 0.05")
    audit.add("Significant comparisons with |Hedges' g| >= 0.8", 379, int((sig & large).sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count q < 0.05 and |g| >= 0.8")
    audit.add("Significant comparisons with |g| < 0.8", 432, int((sig & ~large).sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count q < 0.05 and |g| < 0.8")
    audit.add("Nonsignificant comparisons", 5314, int((~sig).sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count global BH q >= 0.05")
    audit.add("Significant overall (%)", 13.2, 100 * sig.mean(), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "100 * 811 / 6125", .05)
    audit.add("Significant associations meeting large-effect threshold (%)", 46.7, 100 * (sig & large).sum() / sig.sum(), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "100 * 379 / 811", .05)

    changed = int(pd.read_csv(RESULTS / "selection_sensitivity_summary.csv").selected_image_changed.iloc[0])
    alt_sig = alternative.mannwhitney_q_fdr_global.lt(.05)
    primary_set = set(map(tuple, main_results.loc[sig, ["disease", "feature"]].to_numpy()))
    alternative_set = set(map(tuple, alternative.loc[alt_sig, ["disease", "feature"]].to_numpy()))
    shared, primary_only, alternative_only = primary_set & alternative_set, primary_set - alternative_set, alternative_set - primary_set
    union = primary_set | alternative_set
    audit.add("Disease photographs changed", 94, changed, "revision_results/final_20260914/selection_sensitivity_summary.csv", "selected_image_changed")
    audit.add("Selected disease photographs changed (%)", 5.1, 100 * changed / eligible_patients, "revision_results/final_20260914/selection_sensitivity_summary.csv", "100 * changed / 1,850", .05)
    for claim, expected, computed, calculation in [("Alternative-selection significant comparisons", 834, len(alternative_set), "count alternative q < 0.05"),
                                                   ("Shared significant comparisons", 770, len(shared), "set intersection"),
                                                   ("Primary-only significant comparisons", 41, len(primary_only), "primary set difference"),
                                                   ("Alternative-only significant comparisons", 64, len(alternative_only), "alternative set difference"),
                                                   ("Union of significant comparisons", 875, len(union), "set union")]:
        audit.add(claim, expected, computed, "revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz", calculation)
    audit.add("Selection-sensitivity Jaccard overlap (%)", 88.0, 100 * len(shared) / len(union), "revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz", "100 * intersection / union", .05)
    audit.add("Primary significant findings retained (%)", 94.9, 100 * len(shared) / len(primary_set), "revision_results/final_20260914/age_3_9*/significance_testing/stats_each_disease_vs_healthy.csv.gz", "100 * shared / primary significant", .05)

    for layer, label, expected in [("layer1_broad_ethnicity", "Broad", 72), ("layer2_subcategory", "Detailed", 89)]:
        path = RESULTS / f"age_3_9/ethnicity_analysis/{layer}/analysis_A_healthy_kruskal_all_groups.csv.gz"
        frame = pd.read_csv(path)
        audit.add(f"{label} reference omnibus significant", expected, int(truth(frame.significant_fdr).sum()), str(path.relative_to(ROOT)), "count corrected omnibus q < 0.05")
        for analysis, tests, significant in [("B", 36000 if label == "Broad" else 33875, 375 if label == "Broad" else 452),
                                             ("C", 11500 if label == "Broad" else 15000, 479 if label == "Broad" else 633)]:
            filename = "analysis_B_matched_disease_vs_healthy.csv.gz" if analysis == "B" else "analysis_C_all_interactions.csv.gz"
            frame = pd.read_csv(RESULTS / f"age_3_9/ethnicity_analysis/{layer}/{filename}")
            audit.add(f"{label} Analysis {analysis} tests", tests, len(frame), str((RESULTS / f"age_3_9/ethnicity_analysis/{layer}/{filename}").relative_to(ROOT)), "row count")
            audit.add(f"{label} Analysis {analysis} globally significant", significant, int(truth(frame.global_significant_fdr).sum()), str((RESULTS / f"age_3_9/ethnicity_analysis/{layer}/{filename}").relative_to(ROOT)), "count global_significant_fdr=true")

    mapping = pd.read_csv(EVIDENCE / "feature_region_family_mapping.csv")
    enriched = main_results.merge(mapping[["feature", "anatomical_region"]], on="feature", how="left", validate="many_to_one")
    for region, tests, count, percent in [("forehead", 245, 74, 30.2), ("chin / jaw", 245, 66, 26.9), ("global face shape / proportions", 294, 50, 17.0)]:
        rows = enriched[enriched.anatomical_region.eq(region)]
        region_sig = rows.mannwhitney_q_fdr_global.lt(.05)
        audit.add(f"{region} tested comparisons", tests, len(rows), "local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results", "join feature to anatomical region; row count")
        audit.add(f"{region} significant comparisons", count, int(region_sig.sum()), "local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results", "join feature to region; count global q < 0.05")
        audit.add(f"{region} significant (%)", percent, 100 * region_sig.mean(), "local_review/bibm2026_evidence/feature_region_family_mapping.csv + corrected main results", "100 * significant / tested", .05)
    audit.add("Positive significant effects", 490, int((sig & main_results.hedges_g.gt(0)).sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count q < 0.05 and g > 0")
    audit.add("Negative significant effects", 321, int((sig & main_results.hedges_g.lt(0)).sum()), "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", "count q < 0.05 and g < 0")

    recurrence = main_results.assign(significant=sig).groupby("feature").agg(syndromes_tested=("disease", "nunique"), number_significant=("significant", "sum")).reset_index()
    for manuscript_name, feature, count, percent in [("iris_offset_y", "iris_offset_y_l", 24, 49.0), ("forehead_height", "forehead_height", 22, 44.9),
                                                       ("hairline_height", "hairline_height", 22, 44.9), ("philtrum_width", "philtrum_width", 21, 42.9),
                                                       ("philtrum_length", "philtrum_length", 19, 38.8)]:
        row = recurrence[recurrence.feature.eq(feature)].iloc[0]
        audit.add(f"{manuscript_name} recurrence count", count, int(row.number_significant), "corrected main results", "group by exact feature and count q < 0.05", notes="The manuscript shorthand iris_offset_y corresponds to exact feature iris_offset_y_l." if manuscript_name == "iris_offset_y" else "")
        audit.add(f"{manuscript_name} syndromes tested", 49, int(row.syndromes_tested), "corrected main results", "distinct eligible syndrome IDs for exact feature")
        audit.add(f"{manuscript_name} recurrence (%)", percent, 100 * row.number_significant / row.syndromes_tested, "corrected main results", "100 * significant / tested", .05)

    examples = [("Cornelia de Lange syndrome", "philtrum_length", 64, 1.2804, .9891, 1.6064, 7.552e-13),
                ("Noonan syndrome", "inter_pupillary_distance", 34, 1.3189, .8200, 1.8178, 2.881e-6),
                ("Rubinstein-Taybi syndrome", "malar_bulge_mean", 24, -1.0016, -1.4628, -.4923, .002809)]
    for syndrome, feature, n, g, low, high, q in examples:
        row = row_for(main_results, syndrome, feature)
        prefix = f"{syndrome} / {feature} primary"
        for field, expected, computed, tolerance in [("n", n, row.n_disease, 0), ("g", g, row.hedges_g, .00005),
                                                      ("95% CI lower", low, row.hedges_g_ci_95_low, .00005), ("95% CI upper", high, row.hedges_g_ci_95_high, .00005),
                                                      ("global q", q, row.mannwhitney_q_fdr_global, max(abs(q) * .001, 5e-7))]:
            audit.add(f"{prefix} {field}", expected, computed, "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz", f"select exact syndrome and feature; read {field}", tolerance)
    for syndrome, feature, expected_g, expected_low, expected_high, expected_q in [("Cornelia de Lange syndrome", "philtrum_length", 1.3210, None, None, 2.520e-13),
                                                                                  ("Noonan syndrome", "inter_pupillary_distance", 1.32, None, None, None),
                                                                                  ("Rubinstein-Taybi syndrome", "malar_bulge_mean", -.9871, -1.45, -.48, .003719)]:
        row = row_for(alternative, syndrome, feature)
        audit.add(f"{syndrome} / {feature} alternative g", expected_g, row.hedges_g, "revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz", "select exact syndrome and feature; read g", .005 if syndrome != "Cornelia de Lange syndrome" else .00005)
        if expected_low is not None:
            audit.add(f"{syndrome} / {feature} alternative CI lower", expected_low, row.hedges_g_ci_95_low, "revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz", "select exact syndrome and feature; read bootstrap CI", .01)
            audit.add(f"{syndrome} / {feature} alternative CI upper", expected_high, row.hedges_g_ci_95_high, "revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz", "select exact syndrome and feature; read bootstrap CI", .01)
        if expected_q is not None:
            audit.add(f"{syndrome} / {feature} alternative global q", expected_q, row.mannwhitney_q_fdr_global, "revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz", "select exact syndrome and feature; read global q", max(abs(expected_q) * .001, 5e-7))

    broad_b = pd.read_csv(RESULTS / "age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz")
    for group, n, g, q in [("European", 29, 1.1377, 8.343e-4), ("Others", 21, 1.4387, 1.3665e-3)]:
        row = subgroup_row(broad_b, "Cornelia de Lange syndrome", "philtrum_length", group)
        for field, expected, computed, tolerance in [("n", n, row.n_disease, 0), ("g", g, row.hedges_g, .00005), ("global q", q, row.global_fdr_q_value, 5e-7)]:
            audit.add(f"Cornelia de Lange / philtrum_length / {group} {field}", expected, computed, "revision_results/final_20260914/age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_B_matched_disease_vs_healthy.csv.gz", f"select syndrome, feature, and broad group; read {field}", tolerance)

    frame = audit.write()
    bad = frame.status.ne("PASS")
    if bad.any():
        print(frame.loc[bad].to_string(index=False), file=sys.stderr)
        return 1
    print(f"Wrote {OUTPUT_CSV.relative_to(ROOT)} and {OUTPUT_MD.relative_to(ROOT)}: {len(frame)} claims PASS")
    return 0


if __name__ == "__main__": raise SystemExit(main())
