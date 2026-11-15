import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from core import (DEFAULT_OUTPUT, KEYS, ROOT, Checkpoints, Study, baseline, empirical_p, load_baseline,
                  placebo_allocation, permuted_order, rng_for, signature, source_signature, tests_direct, write_json)


def scipy_tests(x, y):
    tests = mannwhitneyu(x, y, axis=0, method='auto', alternative='two-sided')
    return tests.statistic, tests.pvalue


def validate(output):
    study = Study(load_baseline())
    plan = json.loads((output / 'analysis_plan.json').read_text())
    assert source_signature(study.sources) == plan['source_sha256']
    for name, digest in plan['code_sha256'].items():
        assert baseline.sha(ROOT / name) == digest, f'Changed execution code: {name}'
    expected_protected = json.loads((output / 'protected_file_hashes.json').read_text())
    assert baseline.protected_hashes() == expected_protected
    for manifest in [output / 'artifact_sha256.json']:
        if manifest.exists():
            for name, digest in json.loads(manifest.read_text()).items():
                assert baseline.sha(output / name) == digest, f'Changed artifact: {name}'
    primary = study.primary
    reference = pd.read_csv(output / 'reference_resampling_aggregate.csv', float_precision='round_trip')
    permutation = pd.read_csv(output / 'gmdb_permutation_aggregate.csv', float_precision='round_trip')
    placebo = pd.read_csv(output / 'fairface_placebo_summary.csv', float_precision='round_trip')
    states = {stage: Checkpoints(output, stage, signature(plan, stage), True)
              for stage in ['placebo', 'permutation', 'reference-resampling']}
    assert len(reference) == study.total
    assert not reference[KEYS].duplicated().any()
    for name, group in permutation.groupby('null_specification'):
        assert len(group) == study.total and not group[KEYS].duplicated().any()
        p = empirical_p(group.exceedances.to_numpy(), int(group.B.iloc[0]))
        np.testing.assert_allclose(p, group.p_empirical_permutation, rtol=1e-12)
        np.testing.assert_allclose(multipletests(p, method='fdr_bh')[1], group.q_global_BH_empirical, rtol=1e-12)
    reference_state = states['reference-resampling']
    g = np.vstack(reference_state.arrays('g'))
    q = np.vstack(reference_state.arrays('q'))
    np.testing.assert_allclose(np.median(g, axis=0), reference.median_g, rtol=1e-12)
    np.testing.assert_allclose(np.quantile(g, .025, axis=0), reference.g_resampling_q025, rtol=1e-12)
    np.testing.assert_allclose(np.quantile(g, .975, axis=0), reference.g_resampling_q975, rtol=1e-12)
    np.testing.assert_allclose((q < .05).mean(0), reference.significant_frequency)
    np.testing.assert_allclose((np.sign(g) == np.sign(primary.hedges_g.to_numpy())).mean(0), reference.sign_agreement_frequency)
    reps = [0, min(1, len(q) - 1), len(q) - 1]
    for rep in sorted(set(reps)):
        p, direct_g = np.ones(study.total), np.full(study.total, np.nan)
        for i, (race, data) in enumerate(study.races.items()):
            selected = rng_for(plan['seed'], 3, rep, i).choice(len(data['y']), study.sampling['quotas'][race], replace=False)
            for cell_index in data['rows']:
                target = study.cells.disease.iloc[cell_index]
                x, y = data['x'][data['patients'].disease.eq(target).to_numpy()], data['y'][selected]
                where = slice(cell_index * len(study.features), (cell_index + 1) * len(study.features))
                _, p[where] = scipy_tests(x, y)
                direct_g[where] = tests_direct(x, y)[2]
        np.testing.assert_allclose(multipletests(p, method='fdr_bh')[1], q[rep], rtol=1e-10, atol=1e-15)
        np.testing.assert_allclose(direct_g, g[rep], rtol=1e-10, atol=1e-12)
    for rep in sorted(set([0, min(1, len(placebo) - 1), len(placebo) - 1])):
        p, gs = np.ones(study.total), []
        for i, (race, data) in enumerate(study.races.items()):
            _, pairs = placebo_allocation(len(data['y']), data['ns'], study.sampling['quotas'][race], rng_for(plan['seed'], 1, rep, i))
            for cell_index, pair in zip(data['rows'], pairs):
                if pair is None:
                    continue
                x, y = data['y'][pair[0]], data['y'][pair[1]]
                where = slice(cell_index * len(study.features), (cell_index + 1) * len(study.features))
                _, p[where] = scipy_tests(x, y)
                gs.extend(tests_direct(x, y)[2])
        adjusted = multipletests(p, method='fdr_bh')[1]
        assert int((adjusted < .05).sum()) == int(placebo.discoveries.iloc[rep])
        np.testing.assert_allclose(max(np.abs(gs)), placebo.maximum_absolute_g.iloc[rep], rtol=1e-10)
    permutation_state = states['permutation']
    first = permutation_state.state['chunks'][0]
    initial_B = first['stop']
    expected = np.load(permutation_state.directory / first['file'])['counts']
    independently_counted = np.zeros_like(expected)
    for rep in range(initial_B):
        for variant in range(len(expected)):
            for i, (race, data) in enumerate(study.races.items()):
                order = permuted_order(len(data['x']), rng_for(plan['seed'], 2, rep, i, variant), data['blocks'] if variant else None)
                original_codes = np.repeat(np.arange(len(data['labels'])), data['sizes'])
                permuted_codes = np.empty_like(original_codes)
                permuted_codes[order] = original_codes
                sums = np.array([np.bincount(permuted_codes, weights=data['ranks'][:, j], minlength=len(data['labels']))
                                 for j in range(len(study.features))]).T
                n = data['ns'][:, None]
                u = sums[data['positions']] - n * (n + 1) / 2
                center = data['block_center'] if variant else n * (len(data['x']) - n) / 2
                independently_counted[variant, data['features_index']] += (
                    np.abs(u - center) >= np.abs(data['observed_u'] - center) - 1e-10).ravel()
    np.testing.assert_array_equal(expected, independently_counted)
    rng = np.random.default_rng(516)
    null, shifted = [], []
    for _ in range(200):
        y, x = rng.normal(size=(60, 4)), rng.normal(size=(20, 4))
        null.append((multipletests(scipy_tests(x, y)[1], method='fdr_bh')[1] < .05).any())
        shifted.append((multipletests(scipy_tests(x + 2, y)[1], method='fdr_bh')[1] < .05).sum())
    synthetic = {'clearly_fabricated': True, 'seed': 516, 'B': 200, 'features': 4,
                 'null_any_discovery_frequency': float(np.mean(null)), 'shifted_mean_discoveries': float(np.mean(shifted)),
                 'case_n': 20, 'reference_n': 60, 'positive_shift_SD': 2}
    assert synthetic['null_any_discovery_frequency'] < .12 and synthetic['shifted_mean_discoveries'] > 3.5
    write_json(output / 'synthetic_verification.json', synthetic)
    tests = {}
    for name in ['tests.xml', 'regression_tests.xml']:
        if (output / name).exists():
            suites = ET.parse(output / name).getroot()
            tests[name] = {field: sum(int(suite.get(field, 0)) for suite in suites.iter('testsuite'))
                           for field in ['tests', 'errors', 'failures', 'skipped']}
            assert tests[name]['errors'] == tests[name]['failures'] == 0
    summary = {'status': 'PASS', 'planned_family': study.total, 'all_checkpoint_hashes_contiguous': True,
               'input_and_execution_code_hashes_match': True, 'protected_artifacts_unchanged_during_execution': True,
               'direct_SciPy_reference_repetitions': len(set(reps)), 'direct_SciPy_placebo_repetitions': 3,
               'independent_bincount_permutation_repetitions_per_null': initial_B,
               'all_empirical_BH_reconstructed': True, 'all_reference_quantiles_and_frequencies_reconstructed': True,
               'synthetic_verification': synthetic, 'tests': tests}
    write_json(output / 'validation_summary.json', summary)
    path = output / 'RESULTS_AND_LIMITATIONS.md'
    content = path.read_text().split('\n## Execution validation\n')[0]
    path.write_text(content + '\n\n## Execution validation\n\n' + json.dumps(summary, indent=2) + '\n')
    hashes = {str(p.relative_to(output)): baseline.sha(p) for p in sorted(output.rglob('*'))
              if p.is_file() and p.name != 'artifact_sha256.json'}
    write_json(output / 'artifact_sha256.json', hashes)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Independently verify completed private sensitivity outputs')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(validate(baseline.private_output(args.output_dir)), indent=2))
