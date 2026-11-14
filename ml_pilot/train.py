import argparse
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning

from ml_pilot.cv import PRIMARY_SEED, REPEAT_SEEDS, model_spec, nested_cv, permutation_null, tune_model
from ml_pilot.data_access import LABELS, ROOT, S6, PilotError, configuration, load_cohort, sha256, source_plan
from ml_pilot.evaluation import (
    class_tables, coefficient_tables, cohort_summary, confusion_tables, figures, metadata_input, prepare_output,
    stability_table, write_json, write_report, write_table,
)


def parser():
    cli = argparse.ArgumentParser(description="Independent five-cohort PRISM ML pilot; internal local research only.")
    cli.add_argument("--dry-run", action="store_true", help="Audit planned sources without reading patient rows or fitting.")
    cli.add_argument("--cohort", choices=["largest5"], default="largest5")
    cli.add_argument("--authorization", type=Path, help="Private local JSON record of confirmed authorization.")
    cli.add_argument("--seed", type=int, default=PRIMARY_SEED)
    cli.add_argument("--outer-folds", type=int, default=5)
    cli.add_argument("--inner-folds", type=int, default=3)
    cli.add_argument("--no-repeats", action="store_true", help="Omit the two predefined split-stability repetitions.")
    cli.add_argument("--forest", action="store_true", help="Include the declared exploratory random-forest grid.")
    cli.add_argument("--metadata-model", action="store_true", help="Include a separately labeled metadata-only diagnostic.")
    cli.add_argument("--permutations", type=int, default=0, help="0 or at least 100 full nested-CV permutations.")
    cli.add_argument("--permutation-seed", type=int, default=7312026)
    cli.add_argument("--fit-final", action="store_true", help="After evaluation, tune/refit primary logistic on all patients.")
    cli.add_argument("--jobs", type=int, default=1)
    cli.add_argument("--output-dir", type=Path, default=ROOT / "local_review/ml_pilot_20261114")
    return cli


def provenance(config, args, authorization, cohort, models):
    packages = ("scikit-learn", "numpy", "pandas", "scipy", "matplotlib", "joblib", "threadpoolctl")
    sources = {key: sha256(config[key]) for key in ("disease_master_csv", "feature_schema")}
    sources["ACTIVE_ANALYSIS.json"] = sha256(ROOT / "ACTIVE_ANALYSIS.json")
    sources["S6"] = sha256(S6)
    code_files = sorted((ROOT / "ml_pilot").rglob("*.py")) + [ROOT / "ml_pilot/requirements.txt"]
    code = {str(path.relative_to(ROOT)): sha256(path) for path in code_files}
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    return {"status": "RUNNING", "started_at_utc": datetime.now(timezone.utc).isoformat(),
            "study_id": config["analysis_id"], "facekit_commit": config["facekit_commit"],
            "facekit_model_sha256_from_active_config": config["model_sha256"], "source_sha256": sources,
            "code_sha256": code, "git_head": head, "python": platform.python_version(),
            "versions": {name: importlib.metadata.version(name) for name in packages},
            "authorization": authorization, "authorization_record_sha256": sha256(args.authorization),
            "cohort_audit": cohort.audit,
            "protocol": {"cohort": args.cohort, "seed": args.seed, "outer_folds": args.outer_folds,
                         "inner_folds": args.inner_folds, "inner_seed_rule": "outer_seed + 1000 + one_based_outer_fold",
                         "repeat_seeds": [] if args.no_repeats else list(REPEAT_SEEDS),
                         "primary_model": "logistic", "models": models, "grids": {n: model_spec(n, args.seed)[1] for n in models},
                         "selection_metric": "inner balanced accuracy", "primary_metric": "pooled OOF balanced accuracy",
                         "imputation": "training-fold median; all-missing training feature retained as zero",
                         "jobs": args.jobs, "permutations": args.permutations, "permutation_seed": args.permutation_seed,
                         "final_fit_requested": args.fit_final},
            "privacy": {"individual_predictions_persisted": False, "public_release_performed": False,
                        "aggregate_release_requires_pi_review": True, "metadata_minimum_cell_size": 10}}


def verify_sources(manifest, config):
    current = {key: sha256(config[key]) for key in ("disease_master_csv", "feature_schema")}
    current["ACTIVE_ANALYSIS.json"] = sha256(ROOT / "ACTIVE_ANALYSIS.json")
    current["S6"] = sha256(S6)
    if current != manifest["source_sha256"]:
        raise PilotError("BLOCKED_INPUT_CHANGED: source hashes changed during the run; results must not be used.")
    if any(sha256(ROOT / name) != digest for name, digest in manifest["code_sha256"].items()):
        raise PilotError("BLOCKED_CODE_CHANGED: pilot code changed during execution; results must not be used.")


def execute(args, config):
    if args.outer_folds != 5 or args.inner_folds != 3:
        raise PilotError("The real-data primary protocol is fixed to five outer and three inner folds.")
    if args.permutations != 0 and args.permutations < 100:
        raise PilotError("--permutations must be 0 (NOT_RUN) or at least 100.")
    if args.seed in REPEAT_SEEDS and not args.no_repeats:
        raise PilotError("Primary and repeat seeds must be distinct; use --seed 20261114 or --no-repeats.")
    if args.jobs == 0:
        raise PilotError("--jobs must be a nonzero joblib worker count.")
    before = {key: sha256(config[key]) for key in ("disease_master_csv", "feature_schema") if Path(config[key]).is_file()}
    cohort, definitions, authorization = load_cohort(config, args.authorization)
    output = prepare_output(args.output_dir)
    os.environ.setdefault("MPLCONFIGDIR", str(output / "private/matplotlib"))
    models = ["dummy", "logistic", "svm"]
    if args.forest:
        models.append("forest")
    if args.metadata_model:
        models.append("metadata_logistic")
    manifest = provenance(config, args, authorization, cohort, models)
    if any(manifest["source_sha256"][key] != digest for key, digest in before.items()):
        raise PilotError("BLOCKED_INPUT_CHANGED: sources changed while loading; choose a new run.")
    write_json(manifest, output / "run_manifest.json")
    results = []
    try:
        seeds = [args.seed] + ([] if args.no_repeats else list(REPEAT_SEEDS))
        for seed in seeds:
            for name in models:
                X = metadata_input(cohort.metadata) if name == "metadata_logistic" else cohort.X
                result = nested_cv(X, cohort.y, cohort.patient_ids, name, seed, args.outer_folds, args.inner_folds,
                                   args.jobs, progress=lambda message: print(message, flush=True))
                results.append(result)
                print(json.dumps({"completed_model": name, "seed": seed, "pooled_oof": result.score}), flush=True)
        primary = [result for result in results if result.seed == args.seed]
        logistic = next(result for result in primary if result.model == "logistic")
        permutation = {"status": "NOT_RUN", "reason": "--permutations was 0; no p-value calculated"}
        if args.permutations:
            permutation = permutation_null(cohort.X, cohort.y, cohort.patient_ids, logistic.score["balanced_accuracy"],
                                           args.permutations, args.permutation_seed, args.seed, args.outer_folds,
                                           args.inner_folds, args.jobs, progress=lambda message: print(message, flush=True))
        if args.fit_final:
            import joblib
            estimator, selected, inner_score = tune_model(cohort.X, cohort.y, "logistic", args.seed + 2000,
                                                         args.inner_folds, args.jobs)
            destination = output / "private/final_logistic.joblib"
            joblib.dump({"pipeline": estimator, "feature_order": list(cohort.X.columns),
                         "class_mapping": dict(enumerate(LABELS)), "role": "all-data research fit; no independent test"},
                        destination)
            destination.chmod(0o600)
            manifest["final_fit"] = {"status": "FITTED", "selected_parameters": selected,
                                     "inner_selection_score_not_test_performance": inner_score,
                                     "seed": args.seed + 2000, "artifact_sha256": sha256(destination)}
        else:
            manifest["final_fit"] = {"status": "NOT_RUN"}
        verify_sources(manifest, config)
        aggregate = output / "aggregate"
        scores = pd.DataFrame([result.score for result in results])
        write_table(scores, aggregate / "model_scores.csv")
        write_table(pd.concat([result.folds for result in results], ignore_index=True), aggregate / "outer_fold_scores.csv")
        write_table(class_tables(cohort.y, results), aggregate / "class_metrics.csv")
        counts, percent = confusion_tables(cohort.y, logistic.predictions)
        write_table(counts, aggregate / "confusion_matrix_counts.csv")
        write_table(percent, aggregate / "confusion_matrix_percent.csv")
        write_table(cohort_summary(cohort.metadata, cohort.y), aggregate / "cohort_summary.csv")
        write_table(pd.DataFrame({"feature": cohort.X.columns, "missing_cells": cohort.X.isna().sum().to_numpy(),
                                  "n": len(cohort.y)}), aggregate / "feature_missingness.csv")
        features, families = coefficient_tables(logistic, cohort.X.columns.tolist(), definitions)
        write_table(features, aggregate / "coefficient_summary.csv")
        write_table(families, aggregate / "family_coefficient_summary.csv")
        write_table(stability_table(scores), aggregate / "repeat_seed_sensitivity.csv")
        if permutation["status"] == "EXECUTED":
            write_table(pd.DataFrame({"permutation": range(1, args.permutations + 1),
                                      "balanced_accuracy": permutation["null_scores"]}), aggregate / "permutation_null.csv")
        figures(cohort.y, primary, aggregate, permutation)
        write_report(cohort, results, output, manifest, permutation)
        manifest["permutation_check"] = {key: value for key, value in permutation.items() if key != "null_scores"}
        manifest["status"] = "COMPLETED_INTERNAL_VALIDATION"
        manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["aggregate_sha256"] = {path.name: sha256(path) for path in sorted(aggregate.iterdir())}
        write_json(manifest, output / "run_manifest.json")
    except Exception:
        manifest["status"] = "FAILED_DO_NOT_USE_PARTIAL_OUTPUTS"
        write_json(manifest, output / "run_manifest.json")
        raise
    print(json.dumps({"status": manifest["status"], "primary_seed_scores": [result.score for result in primary],
                      "output_directory": str(output.relative_to(ROOT)), "permutation_check": manifest["permutation_check"]},
                     indent=2, allow_nan=False))
    return 0


def main(argv=None):
    os.umask(0o077)
    args = parser().parse_args(argv)
    try:
        config = configuration()
        plan = source_plan(config)
        if args.dry_run:
            plan["authorization"] = "REQUIRED_FOR_PATIENT_LEVEL_EXECUTION"
            print(json.dumps(plan, indent=2))
            return 0 if plan["status"] != "BLOCKED_PRIVATE_DATA_MISSING" else 2
        return execute(args, config)
    except PilotError as error:
        print(str(error), file=sys.stderr)
        return 2
    except ConvergenceWarning:
        print("BLOCKED_MODEL_CONVERGENCE: review numerical convergence before using results.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
