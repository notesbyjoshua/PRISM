import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from research.race_stratified_pitch_sensitivity.design import (
    ARMS, DEFAULT_OUTPUT, METRICS, RACES, allocations, digest, load_study, locked, pose_metrics, protected_hashes, rng_for)
from research.race_stratified_pitch_sensitivity.run_pitch_imbalance import refresh_manifest


def independent_tests(x, y):
    n, m = len(x), len(y)
    ties = np.any(np.diff(np.sort(np.concatenate([x, y]), axis=0), axis=0) == 0, axis=0)
    exact = (~ties) & (min(n, m) <= 8)
    p, u = np.empty(x.shape[1]), np.empty(x.shape[1])
    for take, method in [(exact, 'exact'), (~exact, 'asymptotic')]:
        if take.any():
            result = mannwhitneyu(x[:, take], y[:, take], axis=0, alternative='two-sided', method=method, use_continuity=True)
            p[take], u[take] = result.pvalue, result.statistic
    df = n+m-2
    variance = (np.sum((x-x.mean(0))**2, axis=0)+np.sum((y-y.mean(0))**2, axis=0))/df
    correction = math.exp(math.lgamma(df/2)-.5*math.log(df/2)-math.lgamma((df-1)/2))
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.where(variance > 0, correction*(x.mean(0)-y.mean(0))/np.sqrt(variance), np.nan)
    return u, p, g


def validate(output):
    print('Independently validating sources, selections, numerical checkpoints and all repetition summaries...', flush=True)
    study = load_study()
    plan = json.loads((output/'analysis_plan.json').read_text())
    manifest = json.loads((output/'source_and_execution_manifest.json').read_text())
    assert digest(plan) == manifest['plan_sha256'] and plan == manifest['plan']
    assert locked.source_signature(study['sources']) == plan['source_sha256']
    for field in ['numerical_code_sha256', 'reporting_code_sha256']:
        for name, expected in manifest[field].items():
            assert locked.baseline.sha(ROOT/name) == expected, 'Changed execution/verification code: '+name
    assert protected_hashes(study['sources']) == json.loads((output/'protected_file_hashes.json').read_text())
    for name, expected in manifest['output_sha256'].items():
        assert locked.baseline.sha(output/name) == expected, 'Changed result artifact: '+name
    np.testing.assert_allclose(multipletests(study['primary'].p, method='fdr_bh')[1], study['primary'].q, rtol=1e-12, atol=1e-15)
    saved_support = pd.read_csv(output/'support_audit.csv', float_precision='round_trip')
    pd.testing.assert_frame_equal(saved_support, study['support'], check_dtype=False)
    pd.testing.assert_frame_equal(pd.read_csv(output/'pitch_availability_audit.csv', float_precision='round_trip'), study['availability'], check_dtype=False)
    assert digest(study['supported'].tolist()) == plan['supported_mask_sha256']
    checkpoint = locked.Checkpoints(output, 'paired_placebo', digest(plan), True)
    repetitions = plan['repetitions_per_arm']
    assert checkpoint.state['completed'] == repetitions
    summary = pd.read_csv(output/'repetition_summary.csv', float_precision='round_trip').set_index(['repetition', 'arm'])
    assert summary.index.is_unique and len(summary) == 2*repetitions
    supported_rows = study['support'].loc[study['support'].supported]
    cell_position = {cell: i for i, cell in enumerate(supported_rows.cell_index)}
    where = np.flatnonzero(study['supported'])
    direct_reps = sorted({0, min(1, repetitions-1), repetitions-1})
    for arm in ARMS:
        own = summary.xs(arm, level='arm').discoveries
        direct_reps.append(int(own.idxmax()))
        if own.gt(0).any():
            direct_reps.append(int(own.loc[own.gt(0)].index[0]))
    direct_reps = sorted(set(direct_reps))
    hits = np.zeros((2, 9600), int)
    category_pose = []
    for chunk in checkpoint.state['chunks']:
        with np.load(checkpoint.directory/chunk['file']) as arrays:
            count = chunk['stop']-chunk['start']
            assert arrays['p'].shape == arrays['g'].shape == (count, 2, len(where))
            assert arrays['pose'].shape == (count, 2, len(supported_rows), len(METRICS))
            assert arrays['selection_sha256'].shape == (count, 2)
            assert np.isfinite(arrays['p']).all() and ((arrays['p'] >= 0) & (arrays['p'] <= 1)).all()
            assert np.isfinite(arrays['pose']).all()
            for offset, rep in enumerate(range(chunk['start'], chunk['stop'])):
                for arm_index, arm in enumerate(ARMS):
                    p, g = np.ones(9600), np.full(9600, np.nan)
                    p[where], g[where] = arrays['p'][offset, arm_index], arrays['g'][offset, arm_index]
                    q = multipletests(p, method='fdr_bh')[1]
                    found = (q < .05) & study['supported']
                    hits[arm_index] += found
                    assert (p[~study['supported']] == 1).all() and (q[~study['supported']] == 1).all()
                    saved = summary.loc[(rep, arm)]
                    assert saved.planned_family == 9600 and saved.supported_tests == len(where)
                    assert saved.discoveries == found.sum()
                    assert saved.unique_feature_names_discovered == found.reshape(80, 120).any(0).sum()
                    assert saved.discoveries_abs_g_ge_0_8 == np.sum(found & (np.abs(g) >= .8))
                    np.testing.assert_allclose([saved.minimum_q, saved.minimum_supported_q], [q.min(), q[where].min()], rtol=1e-12, atol=1e-15)
                    selection = []
                    for race_index, race in enumerate(RACES):
                        rows = supported_rows.loc[supported_rows.canonical_race.eq(race)]
                        all_cells = study['support'].canonical_race.eq(race).to_numpy()
                        assert saved[race+'_discoveries'] == found.reshape(80, 120)[all_cells].sum()
                        if rows.empty:
                            continue
                        group = study['groups'][race]
                        _, pairs = allocations(len(group['values']), group['low'], group['high'], rows.n_case.to_numpy(int),
                                                 study['quotas'][race], rng_for(plan['seed'], rep, race_index, arm_index), arm)
                        if rows.random_design.iloc[0] == 'shared_disjoint' and arm == 'random':
                            draw = np.concatenate([pairs[0][1], *[pair[0] for pair in pairs]])
                            assert len(np.unique(draw)) == len(draw)
                        for row, (case, reference) in zip(rows.itertuples(), pairs):
                            assert len(case) == row.n_case and len(reference) == row.n_reference
                            assert not set(case) & set(reference)
                            assert len(set(case)) == len(case) and len(set(reference)) == len(reference)
                            selection.append([row.cell_index, case.tolist(), reference.tolist()])
                            metric = pose_metrics(group['pose'], group['sex'], case, reference)
                            np.testing.assert_allclose(metric, arrays['pose'][offset, arm_index, cell_position[row.cell_index]], rtol=1e-12, atol=1e-12)
                            if arm == 'pitch_imbalanced':
                                assert set(case).issubset(set(group['high'])) and set(reference).issubset(set(group['low']))
                                assert metric[12] > 0 and metric[13] > 0 and metric[15] >= 0 and metric[14] == 0
                            if rep in direct_reps:
                                _, direct_p, direct_g = independent_tests(group['values'][case], group['values'][reference])
                                indices = slice(row.cell_index*120, (row.cell_index+1)*120)
                                np.testing.assert_allclose(p[indices], direct_p, rtol=1e-10, atol=1e-15)
                                np.testing.assert_allclose(g[indices], direct_g, rtol=1e-10, atol=1e-12, equal_nan=True)
                        positions = [cell_position[cell] for cell in rows.cell_index]
                        category_pose.append({'repetition': rep, 'arm': arm, 'canonical_race': race,
                                              **dict(zip(METRICS, arrays['pose'][offset, arm_index, positions].mean(0)))})
                    assert digest(selection) == arrays['selection_sha256'][offset, arm_index]
        print(f'Verified paired repetitions through {chunk["stop"]}/{repetitions}', flush=True)
    expected_pose = pd.DataFrame(category_pose).set_index(['repetition','arm','canonical_race']).sort_index()
    saved_pose = pd.read_csv(output/'pose_gap_by_repetition_category.csv', float_precision='round_trip')
    saved_pose = saved_pose.loc[saved_pose.supported_cells.gt(0)].set_index(['repetition','arm','canonical_race']).sort_index()
    np.testing.assert_allclose(saved_pose[METRICS], expected_pose[METRICS], rtol=1e-12, atol=1e-12)
    recurrence = pd.read_csv(output/'test_discovery_frequencies.csv', float_precision='round_trip')
    for arm_index, arm in enumerate(ARMS):
        own = recurrence.loc[recurrence.arm.eq(arm)]
        assert list(own[locked.KEYS].itertuples(index=False, name=None)) == list(study['primary'][locked.KEYS].itertuples(index=False, name=None))
        np.testing.assert_array_equal(own.discovery_count, hits[arm_index])
        np.testing.assert_allclose(own.loc[own.supported, 'discovery_frequency'], hits[arm_index, where]/repetitions, rtol=1e-12)
        assert own.loc[~own.supported, 'discovery_frequency'].isna().all()
    categories = pd.read_csv(output/'discovery_summary_by_category.csv', float_precision='round_trip')
    overall = pd.read_csv(output/'discovery_summary_overall.csv', float_precision='round_trip')
    for arm in ARMS:
        own = summary.xs(arm, level='arm')
        row = overall.loc[overall.arm.eq(arm)].iloc[0]
        np.testing.assert_allclose([row['mean'],row['median'],row.q25,row.q75,row.q95,row.minimum,row.maximum,row.runs_with_any],
                                  [own.discoveries.mean(),own.discoveries.median(),*own.discoveries.quantile([.25,.75,.95]),own.discoveries.min(),own.discoveries.max(),own.discoveries.gt(0).sum()])
        for race in RACES:
            row = categories.loc[categories.arm.eq(arm) & categories.canonical_race.eq(race)].iloc[0]
            if row.supported_tests:
                np.testing.assert_allclose([row['mean'],row['median'],row.maximum,row.runs_with_any],
                                          [own[race+'_discoveries'].mean(),own[race+'_discoveries'].median(),own[race+'_discoveries'].max(),own[race+'_discoveries'].gt(0).sum()])
            else:
                assert pd.isna(row['mean']) and pd.isna(row.fraction_with_any) and row.status == 'UNSUPPORTED'
    paired = pd.read_csv(output/'paired_discovery_differences.csv')
    a, b = summary.xs('random',level='arm').sort_index(), summary.xs('pitch_imbalanced',level='arm').sort_index()
    np.testing.assert_array_equal(paired.difference_B_minus_A, b.discoveries.to_numpy()-a.discoveries.to_numpy())
    feature_table = pd.read_csv(output/'feature_discovery_frequency.csv')
    for i, arm in enumerate(ARMS):
        np.testing.assert_array_equal(feature_table.loc[feature_table.arm.eq(arm),'total_discovery_tests_across_repetitions'], hits[i].reshape(80,120).sum(0))
    for name in ['pitch_placebo_distribution','pitch_imbalance_diagnostics']:
        with Image.open(output/(name+'.png')) as image:
            assert all(abs(d-300) < .1 for d in image.info['dpi'])
        svg = (output/(name+'.svg')).read_text()
        assert '<image ' not in svg
        assert (output/(name+'.pdf')).read_bytes().startswith(b'%PDF-')
    diagnostic_svg = (output/'pitch_imbalance_diagnostics.svg').read_text()
    assert 'N/A' in diagnostic_svg and '647' in diagnostic_svg
    assert protected_hashes(study['sources']) == json.loads((output/'protected_file_hashes.json').read_text())
    tests = {}
    if (output/'tests.xml').exists():
        suites = ET.parse(output/'tests.xml').getroot()
        tests = {field: sum(int(suite.get(field,0)) for suite in suites.iter('testsuite')) for field in ['tests','failures','errors','skipped']}
        assert tests['errors'] == tests['failures'] == 0
    report = {'status': 'PASS' if repetitions == 1000 else 'PASS_PILOT', 'repetitions_verified_per_arm': repetitions,
              'planned_family': 9600, 'supported_tests': len(where), 'all_2000_BH_families_reconstructed': repetitions == 1000,
              'all_selection_hashes_and_pair_pose_metrics_reproduced': True,
              'independent_SciPy_repetitions_per_arm': direct_reps, 'independent_Hedges_g_math_lgamma': True,
              'unsupported_p_and_q_equal_one': True, 'common_mask_identical': True,
              'all_checkpoint_hashes_and_source_code_verified': True,
              'protected_file_hash_count': manifest['protected_hash_count'], 'all_protected_files_unchanged': True,
              'figures_300_DPI_and_vector': True, 'European_marked_NA': True, 'fabricated_tests': tests,
              'inference_scope': 'New paired FairFace placebo only; prior studies not rerun'}
    locked.write_json(output/'validation_report.json', report)
    status = 'COMPLETED' if repetitions == 1000 else 'PILOT_COMPLETED'
    path = output/'RESULTS_AND_LIMITATIONS.md'
    content = path.read_text().split('\n## Independent validation\n')[0]
    content = content.replace('Status: AWAITING_INDEPENDENT_VALIDATION','Status: '+status).replace('Status: PILOT_AWAITING_VALIDATION','Status: '+status)
    path.write_text(content+'\n## Independent validation\n\n'+json.dumps(report,indent=2)+'\n')
    refresh_manifest(output, status, independently_validated=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Independent source, selection, SciPy and full-family BH verification')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(validate(locked.baseline.private_output(args.output_dir)), indent=2))
