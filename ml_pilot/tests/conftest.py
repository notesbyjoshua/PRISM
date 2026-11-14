import numpy as np
import pandas as pd
import pytest

from ml_pilot.data_access import EXPECTED_COUNTS, FACEKIT_PIN, LABELS


@pytest.fixture
def artificial_cohort():
    features = [f"fixture_feature_{i:03}" for i in range(120)]
    definitions = pd.DataFrame({"feature": features, "order": range(120), "source_commit": FACEKIT_PIN,
                                "feature_family": [f"fixture_family_{i // 4}" for i in range(120)]})
    y = np.repeat(np.arange(5), EXPECTED_COUNTS)
    rng = np.random.default_rng(51)
    X = rng.normal(size=(len(y), 120))
    X[:, :5] += np.eye(5)[y] * 8
    frame = pd.DataFrame(X, columns=features)
    frame["patient_id"] = [f"artificial_subject_{i}" for i in range(len(y))]
    frame["disease"] = np.array(LABELS)[y]
    frame["frontal_ok"] = True
    frame["age_year"] = 5
    frame["age_month"] = 6
    frame["conditional_age_years"] = 5.5
    frame["sex_label"] = "Unknown"
    frame["pose_yaw"] = rng.normal(size=len(y))
    frame["pose_pitch"] = rng.normal(size=len(y))
    frame["pose_roll"] = rng.normal(size=len(y))
    return frame, features, definitions
