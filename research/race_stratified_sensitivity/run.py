import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

from core import (DEFAULT_OUTPUT, HERE, Study, baseline, load_baseline, signature, source_signature, write_json)
from experiments import prepare_fast_paths, run_permutation, run_placebo, run_reference


def main(argv=None):
    parser = argparse.ArgumentParser(description='Three sensitivities for the locked three-race PRISM redesign')
    parser.add_argument('--analysis', choices=['placebo', 'permutation', 'reference-resampling', 'all'], default='all')
    parser.add_argument('--repetitions', type=int, default=1000)
    parser.add_argument('--permutation-repetitions', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=20261115)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    if args.repetitions < 1 or not 1 <= args.permutation_repetitions <= 100000 or args.seed < 0:
        parser.error('Positive repetitions, nonnegative seed, and permutation B<=100000 required')
    output = baseline.private_output(args.output_dir)
    started = time.perf_counter()
    try:
        study = Study(load_baseline())
    except (FileNotFoundError, ValueError, AssertionError, KeyError) as error:
        write_json(output / 'overall_summary.json', {'status': 'BLOCKED', 'reason': str(error)})
        print(f'BLOCKED: {error}', file=sys.stderr)
        return 2
    study.baseline_sample = baseline.sample_reference(study.reference, study.sampling['quotas'], study.sampling['seed'])
    prepare_fast_paths(study)
    plan = study.plan(args.seed, args.repetitions, args.permutation_repetitions)
    protected = baseline.protected_hashes()
    original_plan = output / 'analysis_plan.json'
    if original_plan.exists():
        previous = json.loads(original_plan.read_text())
        if signature(previous, 'all') != signature(plan, 'all'):
            raise ValueError('Existing analysis plan differs; select a fresh output directory')
    write_json(original_plan, plan)
    cohort = {'source_audit': study.audit, 'races': plan['races'], 'all_candidate_cells': len(study.counts),
              'eligible_cells': len(study.cells), 'planned_hypotheses': study.total,
              'same_source_placebo_support': 'Per-race full-disjoint if feasible; per-pair fallback; impossible contrasts p=1',
              'blocked_triggered_cells': sum(item['substantial_imbalance'] for item in study.covariates),
              'primary_and_sensitivity_estimands_differ': True}
    write_json(output / 'cohort_audit.json', cohort)
    study.counts.to_csv(output / 'cohort_candidate_cells.csv', index=False)
    study.cells.to_csv(output / 'locked_eligible_cells.csv', index=False)
    pd.DataFrame(study.covariates).to_csv(output / 'age_sex_exchangeability_audit.csv', index=False)
    mobility = [{'disease': study.cells.disease.iloc[index], 'canonical_race': race,
                 'case_n': int(n), 'movable_case_labels': int(movable), 'immobile_case_labels': int(n - movable)}
                for race, data in study.races.items()
                for index, n, movable in zip(data['rows'], data['ns'], data['movable_case_counts'])]
    pd.DataFrame(mobility).to_csv(output / 'blocked_label_mobility.csv', index=False)
    study.definitions.to_csv(output / 'exact_120_feature_definitions.csv', index=False)
    (output / 'DATA_DICTIONARY.md').write_text((HERE / 'DATA_DICTIONARY.md').read_text())
    preview = {key: plan[key] for key in ['planned_hypotheses', 'eligible_cells', 'races', 'analysis_seeds', 'blocked_enabled']}
    print(json.dumps(preview, indent=2), flush=True)
    if args.dry_run:
        write_json(output / 'dry_run_status.json', {'status': 'PASS', 'sources_available': True, **preview})
        return 0
    commands_path = output / 'execution_commands.json'
    commands = json.loads(commands_path.read_text()) if commands_path.exists() else []
    command = {'argv': vars(args) | {'output_dir': str(output)}, 'started_utc': pd.Timestamp.now(tz='UTC').isoformat()}
    commands.append(command)
    write_json(commands_path, commands)
    if args.analysis in ['placebo', 'all']:
        run_placebo(study, plan, output, args.repetitions, args.resume)
    if args.analysis in ['permutation', 'all']:
        run_permutation(study, plan, output, args.permutation_repetitions, args.resume)
    if args.analysis in ['reference-resampling', 'all']:
        run_reference(study, plan, output, args.repetitions, args.resume)
    if baseline.protected_hashes() != protected:
        raise ValueError('Protected old study/publication artifacts changed during execution')
    if source_signature(study.sources) != plan['source_sha256']:
        raise ValueError('Source or locked primary changed during execution')
    from reporting import report
    summary = report(study, plan, output)
    command['elapsed_seconds'] = time.perf_counter() - started
    command['status'] = summary['status']
    write_json(commands_path, commands)
    write_json(output / 'protected_file_hashes.json', protected)
    hashes = {str(p.relative_to(output)): baseline.sha(p) for p in sorted(output.rglob('*'))
              if p.is_file() and p.name != 'artifact_sha256.json'}
    write_json(output / 'artifact_sha256.json', hashes)
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
