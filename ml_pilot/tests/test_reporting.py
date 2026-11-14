import json
import subprocess

import numpy as np
import pandas as pd
import pytest

from ml_pilot.cv import nested_cv
from ml_pilot.data_access import DISPLAY_LABELS, ROOT, PilotError, validate_cohort
from ml_pilot import evaluation as ev
from ml_pilot import train


def test_five_class_totals_and_order(artificial_cohort):
    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    result = nested_cv(cohort.X, cohort.y, cohort.patient_ids, "dummy")
    counts, percentages = ev.confusion_tables(cohort.y, result.predictions)
    assert counts.true_cohort.tolist() == list(DISPLAY_LABELS)
    assert counts.columns.tolist() == ["true_cohort", *DISPLAY_LABELS]
    assert counts.iloc[:, 1:].to_numpy().sum() == 246
    np.testing.assert_allclose(percentages.iloc[:, 1:].sum(axis=1), 100)
    classes = ev.class_tables(cohort.y, [result])
    assert len(classes) == 5 and classes.support.sum() == 246


def test_exports_reject_identifiers_arrays_and_unapproved_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(ev, "ROOT", tmp_path)
    output = ev.prepare_output(tmp_path / "local_review/test")
    destination = output / "aggregate/model_scores.csv"
    for name in ["patient_id", "image_id", "file", "path", "landmarks", "fixture_feature_000", "prediction"]:
        with pytest.raises(PilotError, match="export rejected"):
            ev.write_table(pd.DataFrame({name: ["artificial"]}), destination)
    with pytest.raises(PilotError, match="nested"):
        ev.write_table(pd.DataFrame({"n": [[1, 2]]}), destination)
    with pytest.raises(PilotError, match="under the ignored"):
        ev.write_table(pd.DataFrame({"n": [246]}), tmp_path / "model_scores.csv")
    ev.write_table(pd.DataFrame({"n": [246]}), destination)
    ev.write_json({"status": "fixture"}, output / "run_manifest.json")
    assert destination.stat().st_mode & 0o777 == 0o600
    assert output.stat().st_mode & 0o777 == 0o700
    assert (output / "private").stat().st_mode & 0o777 == 0o700
    with pytest.raises(PilotError, match="already exists"):
        ev.prepare_output(output)
    link = tmp_path / "local_review/link"
    link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(PilotError, match="Symlinks"):
        ev.prepare_output(link / "escape")


def test_small_metadata_cells_suppress_complements(artificial_cohort):
    frame, features, definitions = artificial_cohort
    frame.loc[0, "sex_label"] = "Female"
    cohort = validate_cohort(frame, features, definitions)
    table = ev.cohort_summary(cohort.metadata, cohort.y)
    sex = table.loc[table.variable.eq("sex_label") & table.cohort.eq(frame.disease.iloc[0])]
    assert len(sex) == 1 and "suppressed" in sex.status.iloc[0]
    assert sex["count"].isna().all()
    assert not {"patient_id", "image_id", "path"} & set(table.columns)


def test_dry_run_never_reads_private_patients_or_fits(monkeypatch, capsys):
    monkeypatch.setattr(train, "configuration", lambda: {})
    monkeypatch.setattr(train, "source_plan", lambda config: {"status": "READY_FOR_AUTHORIZATION"})
    monkeypatch.setattr(train, "load_cohort", lambda *a: pytest.fail("dry run must not read patient data"))
    monkeypatch.setattr(train, "execute", lambda *a: pytest.fail("dry run must not fit"))
    assert train.main(["--dry-run"]) == 0
    assert "balanced_accuracy" not in capsys.readouterr().out


def test_missing_data_cli_returns_blocked_without_metrics(monkeypatch, capsys):
    monkeypatch.setattr(train, "configuration", lambda: {"disease_master_csv": "/absent/master.csv",
                                                        "feature_schema": "/absent/schema.csv"})
    monkeypatch.setattr(train, "source_plan", lambda config: {"status": "BLOCKED_PRIVATE_DATA_MISSING"})
    assert train.main([]) == 2
    result = capsys.readouterr()
    assert "BLOCKED_PRIVATE_DATA_MISSING" in result.err
    assert "balanced_accuracy" not in result.out


def test_complete_fixture_run_writes_only_aggregate_tables(artificial_cohort, tmp_path, monkeypatch):
    from ml_pilot import cv

    frame, features, definitions = artificial_cohort
    cohort = validate_cohort(frame, features, definitions)
    monkeypatch.setattr(ev, "ROOT", tmp_path)
    monkeypatch.setattr(train, "ROOT", tmp_path)
    monkeypatch.setattr(train, "load_cohort", lambda *a: (cohort, definitions, {"fixture_only": True}))
    monkeypatch.setattr(train, "verify_sources", lambda *a: None)
    monkeypatch.setattr(train, "provenance", lambda *a: {
        "study_id": "ARTIFICIAL_FIXTURE_ONLY", "facekit_commit": "not-a-clinical-run",
        "source_sha256": {}, "protocol": {"seed": 20261114, "outer_folds": 5, "inner_folds": 3},
    })
    original = cv.model_spec

    def small_spec(name, seed):
        estimator, grid = original(name, seed)
        return estimator, {key: values[:1] for key, values in grid.items()}

    monkeypatch.setattr(cv, "model_spec", small_spec)
    monkeypatch.setattr(train, "tune_model", cv.tune_model)
    output = tmp_path / "local_review/artificial_fixture"
    args = train.parser().parse_args(["--no-repeats", "--fit-final", "--metadata-model", "--output-dir", str(output)])
    config = {"disease_master_csv": "/absent/artificial.csv", "feature_schema": "/absent/artificial_schema.csv"}
    assert train.execute(args, config) == 0
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["status"] == "COMPLETED_INTERNAL_VALIDATION"
    assert (output / "private/final_logistic.joblib").is_file()
    for path in (output / "aggregate").glob("*.csv"):
        table = pd.read_csv(path)
        assert set(table.columns).issubset(ev.TABLE_COLUMNS[path.name])
        assert not {"patient_id", "image_id", "file", "path", "landmarks", "prediction"} & set(table.columns)
    assert (output / "aggregate/confusion_matrix.png").stat().st_size > 1000
    report = (output / "aggregate/results_report.md").read_text()
    assert "ARTIFICIAL_FIXTURE_ONLY" in report
    assert "NOT_RUN" in report


def test_protected_study_files_have_no_tracked_modifications():
    protected = ["ACTIVE_ANALYSIS.json", "stats_results", "manuscript", "publication", "revision_results"]
    changed = subprocess.check_output(["git", "diff", "--name-only", "HEAD", "--", *protected], cwd=ROOT, text=True)
    assert changed == ""
