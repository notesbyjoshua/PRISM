import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from research.race_stratified_pitch_sensitivity.design import (
    ARMS, DEFAULT_OUTPUT, HERE, SEED, Blocked, analysis_plan, digest, evaluate, load_study, locked, protected_hashes)


def refresh_manifest(output, status, **fields):
    path = output/'source_and_execution_manifest.json'
    manifest = json.loads(path.read_text())
    manifest.update(fields)
    manifest['status'] = status
    manifest['reporting_code_sha256'] = locked.source_signature([HERE/'reporting.py', HERE/'validate.py'])
    manifest['output_sha256'] = {str(p.relative_to(output)): locked.baseline.sha(p) for p in sorted(output.rglob('*'))
                                if p.is_file() and p != path}
    locked.write_json(path, manifest)


def run(args):
    if not 1 <= args.repetitions <= 1000 or args.seed < 0:
        raise ValueError('Use 1–1,000 repetitions and a nonnegative seed; only 1,000 is the final experiment')
    output = locked.baseline.private_output(args.output_dir)
    if (output/'analysis_plan.json').exists() and not args.resume and not args.dry_run:
        raise ValueError('Existing frozen plan requires --resume or a fresh ignored output directory')
    started = time.monotonic()
    print('Validating locked sources, original selection/BH and pose provenance...', flush=True)
    study = load_study()
    plan = analysis_plan(study, args.seed, args.repetitions)
    preflight = {'status': 'PREFLIGHT_PASS', 'planned_family': 9600,
                 'supported_cells': plan['supported_cells'], 'supported_tests': plan['supported_tests'],
                 'repetitions_per_arm': args.repetitions, 'pose_units': 'degrees',
                 'support_by_category': study['support'].groupby('canonical_race').supported.sum().astype(int).to_dict(),
                 'newly_unsupported_after_pitch_split': int(study['support'].newly_unsupported_after_split.sum())}
    if args.dry_run:
        print(json.dumps(preflight, indent=2), flush=True)
        return preflight
    if (output/'analysis_plan.json').exists():
        if json.loads((output/'analysis_plan.json').read_text()) != plan:
            raise ValueError('Frozen plan/source/numerical code changed; do not overwrite, use a fresh private directory')
        protected = json.loads((output/'protected_file_hashes.json').read_text())
        if protected_hashes(study['sources']) != protected:
            raise ValueError('Protected sources or prior results changed; cannot resume')
    else:
        print('Hashing protected sources, prior results, selectors, manuscript and figure releases...', flush=True)
        protected = protected_hashes(study['sources'])
        locked.write_json(output/'protected_file_hashes.json', protected)
        locked.write_json(output/'analysis_plan.json', plan)
        study['support'].to_csv(output/'support_audit.csv', index=False)
        study['availability'].to_csv(output/'pitch_availability_audit.csv', index=False)
        locked.write_json(output/'source_and_execution_manifest.json', {
            'status': 'RUNNING', 'plan': plan, 'plan_sha256': digest(plan),
            'source_sha256': plan['source_sha256'], 'numerical_code_sha256': plan['numerical_code_sha256'],
            'revision': locked.subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'command': sys.argv, 'raw_photographs_read': False, 'GMDB_features_used_in_pseudo_arms': False,
            'protected_hash_count': len(protected), 'output_directory': str(output.relative_to(ROOT))})
    checkpoint = locked.Checkpoints(output, 'paired_placebo', digest(plan), args.resume)
    if checkpoint.state['completed'] > args.repetitions:
        raise ValueError('Checkpoint exceeds the frozen repetition budget')
    ps, gs, poses, selections = [], [], [], []
    for repetition in range(checkpoint.state['completed'], args.repetitions):
        arms = [evaluate(study, args.seed, repetition, arm) for arm in range(len(ARMS))]
        ps.append(np.stack([arm[0] for arm in arms]))
        gs.append(np.stack([arm[1] for arm in arms]))
        poses.append(np.stack([arm[2] for arm in arms]))
        selections.append([arm[3] for arm in arms])
        if (repetition+1) % 25 == 0 or repetition+1 == args.repetitions:
            checkpoint.append(repetition+1, p=np.asarray(ps), g=np.asarray(gs), pose=np.asarray(poses),
                              selection_sha256=np.asarray(selections))
            ps, gs, poses, selections = [], [], [], []
            print(f'Paired placebo: {repetition+1}/{args.repetitions} repetitions per arm', flush=True)
    from research.race_stratified_pitch_sensitivity.reporting import write_reports
    write_reports(study, plan, output, checkpoint)
    if protected_hashes(study['sources']) != protected:
        raise ValueError('Protected files changed during execution')
    status = 'AWAITING_INDEPENDENT_VALIDATION' if args.repetitions == 1000 else 'PILOT_AWAITING_VALIDATION'
    refresh_manifest(output, status, completed_repetitions_per_arm=checkpoint.state['completed'],
                     elapsed_seconds=time.monotonic()-started, protected_sources_unchanged=True)
    result = {**preflight, 'status': status, 'private_output': str(output.relative_to(ROOT))}
    print(json.dumps(result, indent=2), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description='Isolated race-stratified same-source head-pitch imbalance sensitivity')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--repetitions', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=SEED)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    try:
        run(args)
    except (Blocked, FileNotFoundError, ValueError, AssertionError, KeyError) as error:
        audit = {'status': 'BLOCKED', 'reason': str(error), 'outcomes_fabricated': False}
        output = locked.baseline.private_output(args.output_dir)
        if not (output/'analysis_plan.json').exists():
            locked.write_json(output/'BLOCKED_AUDIT.json', audit)
            if getattr(error, 'audit', None) is not None:
                error.audit.to_csv(output/'pitch_availability_audit.csv', index=False)
        print(json.dumps(audit, indent=2), flush=True)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
