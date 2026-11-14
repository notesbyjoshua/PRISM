import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from stats_testing.analysis_config import ROOT, active_configuration

STUDY_ID = "PRISM_120_20261028"
FACEKIT_PIN = "5764bdd4af9360b55a046a16579fb5a890a8182f"
S6 = ROOT / "manuscript/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv"
LABELS = (
    "WILLIAMS-BEUREN SYNDROME; WBS__93e4a782", "Cornelia de Lange syndrome__fb4a389d",
    "Coffin-Siris syndrome__46fda5b4", "Mucopolysaccharidoses__dad2c2d0", "Kabuki syndrome__fa463433",
)
DISPLAY_LABELS = ("Williams-Beuren", "Cornelia de Lange", "Coffin-Siris", "Mucopolysaccharidoses", "Kabuki")
EXPECTED_COUNTS = (67, 64, 43, 37, 35)
METADATA = ("conditional_age_years", "sex_label", "pose_yaw", "pose_pitch", "pose_roll")


class PilotError(ValueError):
    pass


@dataclass
class Cohort:
    X: pd.DataFrame
    y: np.ndarray
    patient_ids: np.ndarray
    metadata: pd.DataFrame
    audit: dict


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configuration():
    try:
        config = active_configuration(require_inputs=False)
    except (ValueError, FileNotFoundError) as error:
        raise PilotError("BLOCKED_ACTIVE_CONFIGURATION: validate the current 120-feature study first.") from error
    if config.get("analysis_id") != STUDY_ID or config.get("facekit_commit") != FACEKIT_PIN:
        raise PilotError("BLOCKED_PROVENANCE_MISMATCH: active study ID or FaceKit pin changed; review the protocol.")
    return config


def source_plan(config):
    sources = {key: {"source": str(Path(config[key]).relative_to(ROOT)), "available": Path(config[key]).is_file()}
               for key in ("disease_master_csv", "feature_schema")}
    available = all(source["available"] for source in sources.values())
    return {"status": "READY_FOR_AUTHORIZATION" if available else "BLOCKED_PRIVATE_DATA_MISSING",
            "study_id": config["analysis_id"], "sources": sources,
            "expected_counts": dict(zip(LABELS, EXPECTED_COUNTS)), "expected_n": sum(EXPECTED_COUNTS),
            "feature_count": 120, "patient_data_read": False, "metrics_computed": False}


def require_authorization(path):
    if path is None or not Path(path).is_file():
        raise PilotError("BLOCKED_AUTHORIZATION_MISSING: supply --authorization with a private JSON record of "
                         "PI/mentor approval covering this classification analysis, GMDB authorization and IRB scope.")
    try:
        record = json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise PilotError("BLOCKED_AUTHORIZATION_INVALID: approval record must be valid JSON.") from error
    if (record.get("analysis") != "PRISM_FIVE_SYNDROME_ML_PILOT" or record.get("approved") is not True
            or record.get("gmdb_research_authorization_confirmed") is not True
            or record.get("irb_scope_confirmed") is not True):
        raise PilotError("BLOCKED_AUTHORIZATION_INVALID: record must explicitly cover the analysis, GMDB and IRB scope.")
    for key in ("approval_reference", "approved_on", "use_conditions"):
        if not isinstance(record.get(key), str) or not record[key].strip():
            raise PilotError("BLOCKED_AUTHORIZATION_INVALID: approval reference, date and use conditions are required.")
    if record.get("internal_local_analysis_only") is not True:
        raise PilotError("BLOCKED_AUTHORIZATION_INVALID: this CLI supports internal local analysis only.")
    return record


def read_csv(path, **kwargs):
    with Path(path).open(newline="") as stream:
        header = next(csv.reader(stream), [])
    if not header or len(header) != len(set(header)):
        raise PilotError("BLOCKED_SCHEMA: empty or duplicated CSV column headers.")
    return pd.read_csv(path, **kwargs)


def validate_schema(features, definitions):
    expected = definitions.sort_values("order").feature.tolist()
    if (len(features) != 120 or len(set(features)) != 120 or features != expected
            or len(definitions) != 120 or definitions.order.tolist() != list(range(120))
            or not definitions.source_commit.eq(FACEKIT_PIN).all()):
        raise PilotError("BLOCKED_SCHEMA: require the exact ordered 120 features and pinned provenance from S6.")
    return features


def validate_cohort(frame, features, definitions, expected_counts=EXPECTED_COUNTS):
    validate_schema(features, definitions)
    required = set(features) | {"patient_id", "disease", "frontal_ok", "conditional_age_years", "age_year", "age_month"}
    if not required.issubset(frame.columns) or frame.columns.duplicated().any():
        raise PilotError("BLOCKED_COHORT_SCHEMA: required feature or cohort-selection columns missing/duplicated.")
    ids = frame.patient_id.astype("string")
    if ids.isna().any() or ids.str.strip().eq("").any() or ids.str.strip().duplicated().any():
        raise PilotError("BLOCKED_PATIENT_IDS: missing, blank or duplicate recorded patient IDs; no IDs logged.")
    if not frame.disease.isin(LABELS).all():
        raise PilotError("BLOCKED_LABELS: cohort contains a missing or unexpected exact label.")
    counts = frame.disease.value_counts().reindex(LABELS, fill_value=0)
    if counts.tolist() != list(expected_counts):
        audit = {label: {"expected": int(n), "observed": int(counts[label])}
                 for label, n in zip(LABELS, expected_counts)}
        raise PilotError("BLOCKED_COHORT_COUNTS: " + json.dumps(audit))
    frontal = frame.frontal_ok.astype("string").str.lower()
    if not frontal.isin(["true", "1", "1.0"]).all():
        raise PilotError("BLOCKED_SELECTION: every selected row must have frontal_ok=true.")
    years = pd.to_numeric(frame.age_year, errors="coerce")
    months = pd.to_numeric(frame.age_month, errors="coerce")
    age = pd.to_numeric(frame.conditional_age_years, errors="coerce")
    valid = years.notna() & months.notna() & years.mod(1).eq(0) & months.mod(1).eq(0) & months.between(0, 11)
    if not (valid.all() and age.ge(3).all() and age.lt(10).all()
            and np.allclose(age, years + months / 12, rtol=0, atol=1e-8)):
        raise PilotError("BLOCKED_SELECTION: original integer year/month and 3 <= age < 10 eligibility must hold.")
    if any(not is_numeric_dtype(frame[f]) or is_bool_dtype(frame[f]) for f in features):
        raise PilotError("BLOCKED_FEATURE_DTYPE: all schema features must have numeric, nonboolean dtypes.")
    X = frame.loc[:, features].astype(float).reset_index(drop=True)
    if np.isinf(X.to_numpy()).any():
        raise PilotError("BLOCKED_NONFINITE_FEATURES: infinite accepted measurements require source review.")
    missing = X.isna().sum()
    audit = {"class_counts": dict(zip(LABELS, map(int, counts))), "n": len(frame), "feature_count": len(features),
             "missing_feature_cells": int(missing.sum()), "features_with_missing_cells": int(missing.gt(0).sum()),
             "all_missing_features": int(missing.eq(len(X)).sum()), "unique_patient_ids": True,
             "selection": "frozen primary, one photo per recorded patient, frontal_ok, ages 3 <= age < 10"}
    y = frame.disease.map(dict(zip(LABELS, range(5)))).to_numpy(dtype=int)
    metadata = frame.reindex(columns=METADATA).reset_index(drop=True)
    return Cohort(X, y, ids.to_numpy(), metadata, audit)


def load_cohort(config, authorization_path):
    if not all(Path(config[key]).is_file() for key in ("disease_master_csv", "feature_schema")):
        raise PilotError("BLOCKED_PRIVATE_DATA_MISSING: restore the configured current private master and feature schema; "
                         "no historical, aggregate, reference or synthetic substitute is permitted.")
    authorization = require_authorization(authorization_path)
    definitions = read_csv(S6)
    schema = read_csv(config["feature_schema"])
    if schema.columns.tolist() != ["feature"]:
        raise PilotError("BLOCKED_SCHEMA: feature schema must have exactly one feature column.")
    features = validate_schema(schema.feature.tolist(), definitions)
    frame = read_csv(config["disease_master_csv"], dtype={"patient_id": "string", "image_id": "string"})
    if not {"patient_id", "disease"}.issubset(frame.columns):
        raise PilotError("BLOCKED_COHORT_SCHEMA: master requires recorded patient ID and disease columns.")
    ids = frame.patient_id
    if (ids.isna().any() or ids.str.strip().eq("").any() or ids.str.strip().duplicated().any()
            or frame.disease.isna().any()):
        raise PilotError("BLOCKED_MASTER_IDENTITIES: the active master has missing/duplicate identities or missing labels.")
    cohort = validate_cohort(frame.loc[frame.disease.isin(LABELS)].copy(), features, definitions)
    return cohort, definitions, authorization
