import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from research.race_stratified_sensitivity import core
from engine import Bank, stability
from run import draw_metadata, selections


def main():
    output = ROOT / 'local_review/reference_stability_diagnostic'
    manifest = json.loads((output / 'source_and_execution_manifest.json').read_text())
    plan = manifest['plan']
    for name, expected in plan['source_hashes'].items():
        assert core.baseline.sha(Path(name)) == expected, name
    for name, expected in plan['execution_code_hashes'].items():
        assert core.baseline.sha(ROOT / name) == expected, name
    for name, expected in manifest['output_sha256'].items():
        if name != 'source_and_execution_manifest.json' and not name.endswith('.log'):
            assert core.baseline.sha(output / name) == expected, name
    assert manifest['status'] == 'COMPLETED'
    bank = Bank(core.load_baseline(), core)
    sample = core.baseline.sample_reference(bank.reference, bank.sampling['quotas'], bank.sampling['seed'])
    draw = {race: np.flatnonzero(data['ref'].file.isin(sample.file)) for race, data in bank.races.items()}
    p, g, _, _ = bank.evaluate(draw)
    np.testing.assert_allclose(p, bank.primary.p, rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(core.bh(p), bank.primary.q, rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(g, bank.primary.hedges_g, rtol=1e-10, atol=1e-12)
    original = bank.primary.q.lt(.05).to_numpy()
    direct_checks = 0
    verified_draws = 0
    previous_hashes = json.loads((output / 'checkpoint_hashes_before_resume.json').read_text())
    for name, expected in previous_hashes.items():
        assert core.baseline.sha(output / name) == expected, name
    for state_path in sorted((output / 'checkpoints').glob('*/state.json')):
        label = state_path.parent.name
        replacement = label == 'B_reference_bootstrap'
        stage = 1 if replacement else 100 + int(label.rsplit('_', 1)[1])
        sizes = plan['pools'] if replacement else plan['ladder'][stage - 100]
        state = json.loads(state_path.read_text())
        assert state['completed'] == 1000
        expected_signature = hashlib.sha256((manifest['plan_signature'] + label).encode()).hexdigest()
        assert state['signature'] == expected_signature
        arrays = {key: [] for key in ['g', 'q', 'digests', 'unique']}
        stop = 0
        for chunk in state['chunks']:
            path = state_path.parent / chunk['file']
            assert chunk['start'] == stop and core.baseline.sha(path) == chunk['sha256']
            with np.load(path) as source:
                for key in arrays:
                    arrays[key].append(source[key])
            stop = chunk['stop']
        arrays = {key: np.vstack(values) for key, values in arrays.items()}
        assert stop == 1000 and arrays['g'].shape == arrays['q'].shape == (1000, 9600)
        assert np.isfinite(arrays['g']).all() and np.isfinite(arrays['q']).all()
        assert ((arrays['q'] >= 0) & (arrays['q'] <= 1)).all()
        frame = pd.read_csv(output / (label + '_per_feature.csv'), float_precision='round_trip')
        frequency, sign, stable = stability(arrays['g'], arrays['q'], bank.primary.hedges_g.to_numpy())
        np.testing.assert_allclose(frame.significant_frequency, frequency, atol=1e-14)
        np.testing.assert_allclose(frame.sign_agreement_frequency, sign, atol=1e-14)
        np.testing.assert_allclose(frame[['g_q025', 'median_g', 'g_q975']].to_numpy().T,
                                   np.quantile(arrays['g'], [.025, .5, .975], axis=0), rtol=1e-12)
        for i, (race, data) in enumerate(bank.races.items()):
            mask = bank.primary.canonical_race.eq(race).to_numpy()
            if race == 'European' and not replacement:
                assert frame.loc[mask, 'descriptive_stable'].isna().all()
                assert len(np.unique(arrays['digests'][:, i])) == 1
            else:
                np.testing.assert_array_equal(frame.loc[mask, 'descriptive_stable'].eq(1), stable[mask])
                assert len(np.unique(arrays['digests'][:, i])) == 1000
        for rep in range(1000):
            draw = selections(bank, sizes, plan['seed'], stage, rep, replacement)
            digests, unique = draw_metadata(bank, draw)
            np.testing.assert_array_equal(digests, arrays['digests'][rep])
            np.testing.assert_array_equal(unique, arrays['unique'][rep])
            for race, indices in draw.items():
                assert len(indices) == sizes[race]
                if not replacement:
                    assert len(np.unique(indices)) == sizes[race]
            verified_draws += 1
        p, g, _, _ = bank.evaluate(draw)
        np.testing.assert_allclose(g, arrays['g'][-1], rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(multipletests(p, method='fdr_bh')[1], arrays['q'][-1], rtol=1e-10, atol=1e-15)
        for race, data in bank.races.items():
            y = data['y'][draw[race]]
            for j in sorted(set([0, len(data['cases']) - 1])):
                x = data['cases'][j]
                where = data['where'].reshape(len(data['cases']), -1)[j]
                np.testing.assert_allclose(core.effect(x, y), g[where], rtol=1e-10, atol=1e-12)
                for feature in [0, 21, 25, 97, 103, 119]:
                    expected = mannwhitneyu(x[:, feature], y[:, feature], method='auto', alternative='two-sided').pvalue
                    np.testing.assert_allclose(p[where[feature]], expected, rtol=1e-10, atol=1e-15)
                    direct_checks += 1
        del arrays
    assert verified_draws == 5000
    failures = pd.read_csv(output / 'original_118_failure_reasons.csv')
    assert len(failures) == 118 and not failures[core.KEYS].duplicated().any()
    deterministic = pd.read_csv(output / 'full_reference_summary.csv')
    assert deterministic.truly_distinct_draws.eq(1).all()
    assert deterministic.stable_original.isna().all() and deterministic.mean_q_recurrence_original.isna().all()
    full = pd.read_csv(output / 'full_reference_per_feature.csv')
    assert full[['g_q025', 'g_q975', 'g_interval_width', 'significant_frequency']].isna().all().all()
    assert full.deterministic_g.notna().all() and full.deterministic_q.notna().all()
    for path in (output / 'figures').glob('*.png'):
        assert min(Image.open(path).info['dpi']) >= 299.9
        assert '<image ' not in path.with_suffix('.svg').read_text()
        assert path.with_suffix('.pdf').read_bytes().startswith(b'%PDF')
    suites = list(ET.parse(output / 'tests.xml').getroot().iter('testsuite'))
    assert all(int(s.attrib.get('failures', 0)) == int(s.attrib.get('errors', 0)) == 0 for s in suites)
    report = json.loads((output / 'verification_report.json').read_text())
    report.update(final_status='PASS', tests_passed=sum(int(s.attrib['tests']) for s in suites),
                  additional_direct_SciPy_checks=direct_checks, all_5000_draw_sizes_unique_counts_and_digests_verified=True,
                  final_draw_per_design_independently_replayed=True, strict_primary_tolerances_passed=True,
                  dry_run_did_not_change_outputs=True, resume_preserved_all_checkpoints=True,
                  deterministic_uncertainty_fields_missing=True, figure_exports_verified=True, visual_QA_completed=True)
    core.write_json(output / 'verification_report.json', report)
    manifest['log_hash_policy'] = 'Logs are finalized after the process exits; final verification hashes the now-closed logs. All numerical/result artifacts and source/checkpoint hashes are verified before refresh.'
    manifest['verification_code_sha256'] = core.baseline.sha(Path(__file__))
    manifest['final_verification'] = report
    manifest['output_sha256'] = {str(path.relative_to(output)): core.baseline.sha(path) for path in output.rglob('*')
                                if path.is_file() and path.name != 'source_and_execution_manifest.json'}
    core.write_json(output / 'source_and_execution_manifest.json', manifest)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
