#!/usr/bin/env python3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT = ROOT / "publication/CANDIDATE_PHENOTYPE_EVIDENCE.csv"


def main():
    candidates = pd.read_csv(EVIDENCE / "main_candidate_review.csv")
    mapping = pd.read_csv(EVIDENCE / "feature_region_family_mapping.csv")
    demographic = pd.read_csv(EVIDENCE / "candidate_demographic_support.csv")
    interactions = pd.read_csv(EVIDENCE / "candidate_interaction_support.csv")

    demographic["support"] = demographic.apply(lambda r: f"{r.demographic_layer}/{r.group}: n={r.n_disease}, g={r.hedges_g:.4g}, q={r.global_fdr_q_value:.4g}, {r.evidence_level}", axis=1)
    demographic = demographic.groupby(["syndrome", "feature"], as_index=False).support.agg("; ".join).rename(columns={"support": "demographic_subgroup_support"})
    interactions["support"] = interactions.apply(lambda r: f"{r.demographic_layer}/{r.group_1} vs {r.group_2}: beta={r.interaction_beta_group2_minus_group1:.4g}, "
                                                            f"95% CI {r.interaction_ci_lower:.4g} to {r.interaction_ci_upper:.4g}, q={r.global_fdr_q_value:.4g}", axis=1)
    interactions = interactions.groupby(["syndrome", "feature"], as_index=False).support.agg("; ".join).rename(columns={"support": "formal_interaction_evidence"})
    columns = ["feature", "measurement_description", "normalization_or_formula_context"]
    frame = candidates.merge(mapping[columns], on="feature", validate="many_to_one")
    frame = frame.merge(demographic, on=["syndrome", "feature"], how="left").merge(interactions, on=["syndrome", "feature"], how="left")
    frame["absolute_delta_g"] = (frame.hedges_g_alternative - frame.hedges_g).abs()
    frame["literature_classification"] = "unreviewed"
    frame["clinical_review_status"] = "HUMAN_INPUT_REQUIRED - literature and clinical review"
    frame["demographic_subgroup_support"] = frame.demographic_subgroup_support.fillna("none in current evidence shortlist")
    frame["formal_interaction_evidence"] = frame.formal_interaction_evidence.fillna("none in current evidence shortlist")
    frame = frame.rename(columns={"hedges_g": "primary_g", "hedges_g_ci_95_low": "primary_ci_95_lower", "hedges_g_ci_95_high": "primary_ci_95_upper",
                                  "global_q": "primary_global_q", "hedges_g_alternative": "alternative_g", "alternative_global_q": "alternative_q",
                                  "measurement_description": "anatomical_definition", "normalization_or_formula_context": "normalization_reference"})
    frame[["syndrome", "feature", "n_disease", "primary_g", "primary_ci_95_lower", "primary_ci_95_upper", "primary_global_q", "alternative_g",
           "alternative_q", "absolute_delta_g", "demographic_subgroup_support", "formal_interaction_evidence", "anatomical_definition", "normalization_reference",
           "literature_classification", "clinical_review_status"]].to_csv(OUTPUT, index=False)
    print(f"Generated {len(frame)} candidate evidence rows; all literature and clinical classifications remain unresolved")


if __name__ == "__main__": main()
