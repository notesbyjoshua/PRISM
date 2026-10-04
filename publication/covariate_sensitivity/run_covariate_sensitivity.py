#!/usr/bin/env python3
import ast
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "publication/covariate_sensitivity"
DISEASE_PATH = ROOT / "data/master_dataset_disease.csv"
REFERENCE_PATH = ROOT / "data/master_dataset_healthy.csv"
MAPPING_PATH = ROOT / "stats_testing/ancestry_analysis.py"
CANDIDATES = [
    ("Cornelia de Lange syndrome", "philtrum_length"),
    ("Noonan syndrome", "inter_pupillary_distance"),
    ("Rubinstein-Taybi syndrome", "malar_bulge_mean"),
]
MAPPING_NAMES = {"norm_label", "map_gmdb_layer1", "map_gmdb_layer2", "map_fairface_layer1", "map_fairface_layer2"}


def mapping_functions():
    tree = ast.parse(MAPPING_PATH.read_text())
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in MAPPING_NAMES]
    if {node.name for node in nodes} != MAPPING_NAMES: raise ValueError("Required historical demographic mappings are unavailable")
    scope = {"pd": pd, "np": np, "re": __import__("re")}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(MAPPING_PATH), "exec"), scope)
    return scope


def truth(values):
    return values.astype(str).str.strip().str.lower().eq("true")


def normalized_sex(values):
    mapped = values.astype("string").str.strip().str.lower().map({"male": "Male", "m": "Male", "female": "Female", "f": "Female"})
    return mapped.astype("string")


def disease_age(frame):
    year = pd.to_numeric(frame.age_year, errors="coerce")
    month = pd.to_numeric(frame.age_month, errors="coerce")
    valid = year.notna() & month.notna() & year.ge(0) & month.between(0, 11) & year.mod(1).eq(0) & month.mod(1).eq(0)
    return (year + month / 12).where(valid & ~(year.eq(0) & month.eq(0)))


def select_cohorts():
    disease = pd.read_csv(DISEASE_PATH, dtype={"patient_id": "string", "image_id": "string"}, low_memory=False)
    reference = pd.read_csv(REFERENCE_PATH, low_memory=False)
    ages = disease_age(disease).rename("conditional_age_years")
    pose = disease[["pose_yaw", "pose_pitch", "pose_roll"]].apply(pd.to_numeric, errors="coerce").abs()
    pose_score = (pose / [15, 15, 10]).max(axis=1).where(np.isfinite(pose).all(axis=1), np.inf).rename("pose_score")
    disease = pd.concat([disease, ages, pose_score], axis=1).copy()
    eligible = truth(disease.frontal_ok) & disease.patient_id.notna() & disease.conditional_age_years.ge(3) & disease.conditional_age_years.lt(10)
    disease = disease.loc[eligible].sort_values(["patient_id", "pose_score", "image_id"], kind="stable").drop_duplicates("patient_id")
    reference_key = reference.file.astype("string").str.strip().str.replace("\\", "/", regex=False).rename("reference_key")
    reference = pd.concat([reference, reference_key], axis=1).copy()
    reference = reference.loc[truth(reference.frontal_ok)].drop_duplicates("reference_key", keep="first")
    return disease.copy(), reference.copy()


def add_harmonized_covariates(disease, reference):
    maps = mapping_functions()
    disease["sex_label"] = normalized_sex(disease.gender)
    reference["sex_label"] = normalized_sex(reference.gender)
    disease["broad_group"] = disease.image_dataset_ethnicity_category.map(maps["map_gmdb_layer1"])
    reference["broad_group"] = reference.race.map(maps["map_fairface_layer1"])
    disease["detailed_group"] = [maps["map_gmdb_layer2"](a, b) for a, b in zip(
        disease.image_dataset_ethnicity_category, disease.image_dataset_ethnicity_sub_category)]
    reference["detailed_group"] = reference.race.map(maps["map_fairface_layer2"])
    for frame in (disease, reference):
        for column in ["pose_yaw", "pose_pitch", "pose_roll"]: frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return disease, reference


def availability_row(label, disease_field, reference_field, disease_values, reference_values, comparable, usable, reason):
    disease_missing, reference_missing = int(disease_values.isna().sum()), int(reference_values.isna().sum())
    return {
        "covariate": label, "disease_field": disease_field, "reference_field": reference_field,
        "disease_n": len(disease_values), "disease_missing_n": disease_missing,
        "disease_missing_fraction": disease_missing / len(disease_values),
        "disease_usable_value_or_category_count": int(disease_values.nunique(dropna=True)),
        "reference_n": len(reference_values), "reference_missing_n": reference_missing,
        "reference_missing_fraction": reference_missing / len(reference_values),
        "reference_usable_value_or_category_count": int(reference_values.nunique(dropna=True)),
        "semantics_comparable": comparable, "usable_in_sensitivity_model": usable, "decision_reason": reason,
    }


def availability(disease, reference):
    rows = []
    rows.append(availability_row("stored frontal-quality information", "frontal_ok", "frontal_ok", disease.frontal_ok, reference.frontal_ok,
                                 "yes as a common eligibility rule", "no as a covariate",
                                 "Both analytic cohorts require a passing stored flag, leaving no modeled variation."))
    for column in ["pose_yaw", "pose_pitch", "pose_roll"]:
        rows.append(availability_row(column.removeprefix("pose_"), column, column, disease[column], reference[column], "yes", "yes",
                                     "Numeric degrees from the same FaceKit pose fields; modeled continuously."))
    rows.append(availability_row("sex / gender label", "gender", "gender", disease.sex_label, reference.sex_label,
                                 "limited; harmonized dataset labels", "yes",
                                 "Only common Male/Female labels are used; labels are dataset annotations."))
    rows.append(availability_row("broad demographic group", "image_dataset_ethnicity_category", "race", disease.broad_group, reference.broad_group,
                                 "limited; established publication mapping", "yes",
                                 "Uses the existing broad harmonization and excludes unmapped labels."))
    rows.append(availability_row("detailed demographic group", "image_dataset_ethnicity_category + image_dataset_ethnicity_sub_category", "race",
                                 disease.detailed_group, reference.detailed_group, "limited; established publication mapping", "no",
                                 "Comparable mapped labels exist, but sparse syndrome-specific cells make the broad layer more stable for three models."))
    rows.append(availability_row("age information / age category", "age_year + age_month", "age", disease.conditional_age_years,
                                 reference.age.astype("string"), "no for exact age", "no",
                                 "Disease age can be conditionally reconstructed; FairFace supplies the 3-9 category, not exact age."))
    return pd.DataFrame(rows)


def continuous_rows(variable, disease_values, reference_values):
    rows = []
    for cohort, values in [("disease", disease_values), ("reference", reference_values)]:
        values = pd.to_numeric(values, errors="coerce")
        finite = values[np.isfinite(values)]
        rows.append({
            "variable": variable, "summary_type": "continuous", "level_or_statistic": "distribution",
            "cohort": cohort, "n": len(values), "missing_n": int(values.isna().sum()), "count": len(finite),
            "percent": np.nan, "mean": finite.mean(), "standard_deviation": finite.std(ddof=1), "median": finite.median(),
            "q1": finite.quantile(.25), "q3": finite.quantile(.75), "minimum": finite.min(), "maximum": finite.max(),
        })
    return rows


def categorical_rows(variable, disease_values, reference_values):
    rows = []
    levels = sorted(set(disease_values.dropna().astype(str)) | set(reference_values.dropna().astype(str))) + ["<missing>"]
    for cohort, values in [("disease", disease_values), ("reference", reference_values)]:
        displayed = values.astype("string").fillna("<missing>")
        for level in levels:
            count = int(displayed.eq(level).sum())
            rows.append({
                "variable": variable, "summary_type": "categorical", "level_or_statistic": level, "cohort": cohort,
                "n": len(values), "missing_n": int(values.isna().sum()), "count": count, "percent": count / len(values),
                "mean": np.nan, "standard_deviation": np.nan, "median": np.nan, "q1": np.nan, "q3": np.nan,
                "minimum": np.nan, "maximum": np.nan,
            })
    return rows


def balance(disease, reference):
    rows = []
    for column in ["pose_yaw", "pose_pitch", "pose_roll"]: rows += continuous_rows(column, disease[column], reference[column])
    rows += categorical_rows("sex_label", disease.sex_label, reference.sex_label)
    rows += categorical_rows("broad_demographic_group", disease.broad_group, reference.broad_group)
    rows += categorical_rows("disease_exact_age_vs_reference_age_category", pd.cut(
        disease.conditional_age_years, [3, 10], right=False, labels=["3-9"] ).astype("string"), reference.age.astype("string"))
    return pd.DataFrame(rows)


def clean_disease_name(values):
    return values.astype(str).str.replace(r"__[0-9a-f]{8}$", "", regex=True)


def adjusted_models(disease, reference):
    disease = disease.copy()
    disease["clean_disease"] = clean_disease_name(disease.disease)
    rows = []
    for syndrome, feature in CANDIDATES:
        case = disease.loc[disease.clean_disease.eq(syndrome)].copy()
        control = reference.copy()
        case["disease_indicator"], control["disease_indicator"] = 1, 0
        combined = pd.concat([case, control], ignore_index=True, sort=False)
        combined["value"] = pd.to_numeric(combined[feature], errors="coerce")
        required = ["value", "disease_indicator", "sex_label", "broad_group", "pose_yaw", "pose_pitch", "pose_roll"]
        model_data = combined.dropna(subset=required).copy()
        model_data = model_data[model_data.sex_label.isin(["Female", "Male"])]
        model_data = model_data[model_data.broad_group.isin(["European", "African", "Asian", "Others"])]
        if model_data.disease_indicator.value_counts().min() < 3: raise ValueError(f"Too few complete cases for {syndrome} / {feature}")
        model = smf.ols("value ~ disease_indicator + C(sex_label) + C(broad_group) + pose_yaw + pose_pitch + pose_roll", data=model_data).fit(cov_type="HC3")
        beta, se, p = model.params["disease_indicator"], model.bse["disease_indicator"], model.pvalues["disease_indicator"]
        ci = model.conf_int().loc["disease_indicator"]
        rows.append({
            "predefined_candidate_order": len(rows) + 1, "syndrome": syndrome, "feature": feature,
            "model": "OLS with HC3 robust standard errors", "estimand": "adjusted disease-vs-reference mean difference in the historical normalized feature",
            "included_observations": "one selected age-eligible disease image per patient plus deduplicated 3-9 FairFace reference images with complete modeled covariates",
            "missing_covariate_handling": "complete-case analysis", "covariates": "sex label; broad demographic group; yaw; pitch; roll",
            "n_total": int(model.nobs), "n_disease": int(model_data.disease_indicator.sum()),
            "n_reference": int((1 - model_data.disease_indicator).sum()), "adjusted_disease_coefficient": beta,
            "hc3_standard_error": se, "ci_95_lower": ci.iloc[0], "ci_95_upper": ci.iloc[1], "raw_p": p,
            "model_rank": int(model.model.rank), "design_columns": int(model.model.exog.shape[1]),
            "condition_number": float(model.condition_number), "exploratory_family": "three predefined highlighted candidate models",
        })
    frame = pd.DataFrame(rows)
    frame["bh_q_across_three_predefined_models"] = multipletests(frame.raw_p, method="fdr_bh")[1]
    return frame


def feasibility_report(frame, disease, reference):
    usable = frame.loc[frame.usable_in_sensitivity_model.eq("yes"), "covariate"].tolist()
    lines = [
        "# Covariate feasibility", "",
        "This audit evaluates whether fields recorded in both source datasets can support a narrow sensitivity analysis.",
        "It does not treat harmonized dataset labels as genetic ancestry and does not treat the FairFace `3-9` category as",
        "an exact age.", "", "## Analytic cohorts", "",
        f"- Disease: {len(disease):,} selected patients aged 3 to under 10 under the frozen selection rule.",
        f"- Reference: {len(reference):,} deduplicated FairFace images with a passing stored frontal-quality flag.",
        "", "## Decision", "",
        "A limited sensitivity analysis is feasible for the three predefined highlighted candidates. It uses " +
        ", ".join(usable) + ". Stored frontal quality is an inclusion rule rather than a covariate because all retained",
        "records pass it. Detailed demographic labels are not used because syndrome-specific cells are sparse.", "",
        "Exact-age adjustment is not feasible. Disease ages can be reconstructed conditionally from year and month fields,",
        "whereas the reference dataset supplies the categorical label `3-9`. The sensitivity models therefore retain the",
        "frozen childhood cohort but do not include age as a numeric or matched covariate.", "", "## Field decisions", "",
        "| Covariate | Comparable semantics | Use in model | Reason |", "|---|---|---|---|",
    ]
    for row in frame.itertuples(): lines.append(f"| {row.covariate} | {row.semantics_comparable} | {row.usable_in_sensitivity_model} | {row.decision_reason} |")
    (OUTPUT / "COVARIATE_FEASIBILITY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def balance_report(balance_frame):
    lines = [
        "# Cohort balance report", "",
        "This report describes recorded covariates in the frozen childhood disease cohort and the deduplicated FairFace",
        "reference cohort. It does not report a single balance score because disease status and source dataset are perfectly",
        "confounded, and exact age is not comparable across sources.", "", "## Pose", "",
        "| Variable | Cohort | N available | Mean | SD | Median | Q1 | Q3 |", "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in balance_frame.loc[balance_frame.summary_type.eq("continuous")].itertuples():
        lines.append(f"| {row.variable} | {row.cohort} | {int(row['count']) if isinstance(row, pd.Series) else int(row.count)} | {row.mean:.3f} | "
                     f"{row.standard_deviation:.3f} | {row.median:.3f} | {row.q1:.3f} | {row.q3:.3f} |")
    lines += ["", "## Categorical distributions", "",
              "Sex-label, broad demographic-label, and available age-category counts are recorded in `cohort_balance.csv`.",
              "These are dataset annotations. The disease exact-age field and reference age category are shown separately",
              "and are not presented as directly balanced measurements."]
    (OUTPUT / "COHORT_BALANCE_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def sensitivity_report(models):
    lines = [
        "# Covariate sensitivity report", "", "## Predefined analysis", "",
        "Before fitting any adjusted model, the analysis set was fixed to three manuscript-highlighted examples:", "",
    ]
    lines += [f"- {syndrome} / `{feature}`" for syndrome, feature in CANDIDATES]
    lines += [
        "", "Each exploratory model is ordinary least squares with HC3 robust standard errors. The estimand is the adjusted",
        "disease-versus-reference mean difference in the historical normalized feature. Models include one selected image",
        "per age-eligible disease patient and deduplicated reference images with complete sex-label, broad demographic-label,",
        "yaw, pitch, and roll fields. Missing modeled covariates use complete-case exclusion. The three disease coefficients",
        "form one predefined exploratory family and receive Benjamini–Hochberg correction across three tests.", "",
        "The original Mann–Whitney analyses remain primary. These models neither replace them nor change their correction family.",
        "", "## Results", "", "| Syndrome | Feature | Disease N | Reference N | Adjusted coefficient | HC3 95% CI | Raw p | Exploratory q |",
        "|---|---|---:|---:|---:|---|---:|---:|",
    ]
    for row in models.itertuples():
        lines.append(f"| {row.syndrome} | `{row.feature}` | {row.n_disease} | {row.n_reference} | {row.adjusted_disease_coefficient:.6g} | "
                     f"[{row.ci_95_lower:.6g}, {row.ci_95_upper:.6g}] | {row.raw_p:.4g} | {row.bh_q_across_three_predefined_models:.4g} |")
    lines += [
        "", "These analyses assess sensitivity to recorded covariates but cannot separate disease effects from unmeasured",
        "acquisition, curation, or source-dataset effects because disease status and dataset source remain perfectly confounded.",
        "", "The adjusted coefficients should be read only as a recorded-covariate sensitivity check. They do not establish",
        "replication, eliminate dataset effects, or support new clinical phenotype claims.",
    ]
    (OUTPUT / "COVARIATE_SENSITIVITY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    disease, reference = add_harmonized_covariates(*select_cohorts())
    availability_frame = availability(disease, reference)
    balance_frame = balance(disease, reference)
    models = adjusted_models(disease, reference)
    availability_frame.to_csv(OUTPUT / "covariate_availability.csv", index=False)
    balance_frame.to_csv(OUTPUT / "cohort_balance.csv", index=False)
    models.to_csv(OUTPUT / "candidate_adjusted_models.csv", index=False)
    feasibility_report(availability_frame, disease, reference)
    balance_report(balance_frame)
    sensitivity_report(models)
    print(f"Wrote covariate feasibility, balance, and {len(models)} predefined adjusted models")


if __name__ == "__main__": main()
