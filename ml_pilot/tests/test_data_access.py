import json

import numpy as np
import pytest

from ml_pilot.data_access import (
    EXPECTED_COUNTS, LABELS, PilotError, load_cohort, require_authorization, validate_cohort, validate_schema,
)


def test_exact_schema_order_and_no_metadata(artificial_cohort):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    assert cohort.X.columns.tolist() == features
    assert cohort.X.shape == (246, 120)
    assert np.bincount(cohort.y).tolist() == list(EXPECTED_COUNTS)
    assert cohort.audit["class_counts"] == dict(zip(LABELS, EXPECTED_COUNTS))
    for bad in [features[::-1], features[:-1], features[:-1] + ["age_year"], features + ["patient_id"]]:
        with pytest.raises(PilotError, match="BLOCKED_SCHEMA"):
            validate_schema(bad, definitions)


@pytest.mark.parametrize("change,match", [
    ("duplicate", "PATIENT_IDS"), ("missing_id", "PATIENT_IDS"), ("blank_id", "PATIENT_IDS"),
    ("unknown_label", "LABELS"), ("missing_label", "LABELS"), ("count", "COUNTS"),
    ("frontal", "SELECTION"), ("age", "SELECTION"), ("month", "SELECTION"),
    ("string_feature", "DTYPE"), ("inf", "NONFINITE"),
])
def test_rejects_invalid_cohort(artificial_cohort, change, match):
    frame, features, definitions = artificial_cohort
    if change == "duplicate":
        frame.loc[1, "patient_id"] = frame.loc[0, "patient_id"]
    elif change in ["missing_id", "blank_id"]:
        frame.loc[0, "patient_id"] = None if change == "missing_id" else " "
    elif change in ["unknown_label", "missing_label"]:
        frame.loc[0, "disease"] = "unrecognized" if change == "unknown_label" else None
    elif change == "count":
        frame = frame.iloc[:-1]
    elif change == "frontal":
        frame.loc[0, "frontal_ok"] = False
    elif change == "age":
        frame.loc[0, "conditional_age_years"] = 10
    elif change == "month":
        frame.loc[0, "age_month"] = 12
    elif change == "string_feature":
        frame[features[0]] = "0"
    else:
        frame.loc[0, features[0]] = np.inf
    with pytest.raises(PilotError, match=match):
        validate_cohort(frame, features, definitions)


def test_missing_cells_retained_for_fold_imputation(artificial_cohort):
    frame, features, definitions = artificial_cohort
    frame.loc[0, features[0]] = np.nan
    result = validate_cohort(frame, features, definitions)
    assert result.audit["missing_feature_cells"] == 1
    assert np.isnan(result.X.iloc[0, 0])


def test_missing_private_data_never_falls_back(tmp_path):
    config = {"disease_master_csv": str(tmp_path / "absent.csv"), "feature_schema": str(tmp_path / "absent_schema.csv")}
    with pytest.raises(PilotError, match="BLOCKED_PRIVATE_DATA_MISSING"):
        load_cohort(config, None)


def test_authorization_required_before_private_read(tmp_path, monkeypatch):
    source = tmp_path / "not_patient_data.csv"
    source.write_text("fixture\n")
    config = {"disease_master_csv": str(source), "feature_schema": str(source)}
    monkeypatch.setattr("ml_pilot.data_access.read_csv", lambda *a, **k: pytest.fail("read before authorization"))
    with pytest.raises(PilotError, match="BLOCKED_AUTHORIZATION_MISSING"):
        load_cohort(config, None)
    approval = tmp_path / "approval.json"
    approval.write_text(json.dumps({"approved": True}))
    with pytest.raises(PilotError, match="BLOCKED_AUTHORIZATION_INVALID"):
        require_authorization(approval)
