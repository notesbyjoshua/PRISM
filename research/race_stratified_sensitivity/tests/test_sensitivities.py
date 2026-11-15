import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.special import gamma
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import core
from experiments import (ReferenceBank, permutation_counts, prepare_fast_paths, run_permutation,
                         run_placebo, run_reference, support_table)


def synthetic_study():
    rng = np.random.default_rng(451)
    features = [f'SYNTHETIC_feature_{i}' for i in range(4)]
    patients, controls = [], []
    for race in core.RACES:
        for syndrome in ['SYNTHETIC_A', 'SYNTHETIC_B']:
            for i in range(6):
                patients.append({'patient_id': f'SYNTHETIC_{race}_{syndrome}_{i}',
                                 'image_id': f'SYNTHETIC_image_{race}_{syndrome}_{i}', 'disease': syndrome,
                                 'canonical_race': race, 'image_dataset_ethnicity_category': race,
                                 'sex_label': 'Male' if i % 2 else 'Female', 'conditional_age_years': 3.5 + i,
                                 **dict(zip(features, rng.normal(size=4)))})
        for i in range(40):
            controls.append({'file': f'SYNTHETIC/{race}/{i:03d}', 'canonical_race': race,
                             **dict(zip(features, rng.normal(size=4)))})
    disease, reference = pd.DataFrame(patients), pd.DataFrame(controls)
    counts = core.baseline.cohort_counts(disease)
    sampled = core.baseline.sample_reference(reference, dict.fromkeys(core.RACES, 20), 1)
    primary = []
    for cell in counts.sort_values(['disease', 'canonical_race']).itertuples():
        if not cell.eligible:
            continue
        x = disease.loc[disease.disease.eq(cell.disease) & disease.canonical_race.eq(cell.canonical_race), features]
        y = sampled.loc[sampled.canonical_race.eq(cell.canonical_race), features]
        u, p, g = core.tests_direct(x, y)
        for j, feature in enumerate(features):
            primary.append({'disease': cell.disease, 'canonical_race': cell.canonical_race, 'feature': feature,
                            'case_n': len(x), 'reference_n': len(y), 'mannwhitney_u': u[j], 'p': p[j],
                            'hedges_g': g[j], 'supported': True})
    primary = pd.DataFrame(primary)
    primary['q'] = core.bh(primary.p)
    sampling = {'quotas': dict.fromkeys(core.RACES, 20), 'X': 60, 'seed': 1}
    study = core.Study((disease, reference, features, pd.DataFrame(), {}, counts, primary, sampling, []))
    study.baseline_sample = sampled
    prepare_fast_paths(study)
    return study


@pytest.mark.parametrize('capacity,sizes,reference_n,design,supported', [
    (50, [5, 6, 8], 20, 'shared_disjoint', 3),
    (30, [8, 9], 20, 'per_pair_disjoint_fallback', 2),
    (20, [5, 6], 20, 'per_pair_disjoint_fallback', 0),
])
def test_disjoint_allocation_capacity_counts_and_no_replacement(capacity, sizes, reference_n, design, supported):
    actual, pairs = core.placebo_allocation(capacity, sizes, reference_n, np.random.default_rng(17))
    assert actual == design and sum(pair is not None for pair in pairs) == supported
    for size, pair in zip(sizes, pairs):
        if pair is None:
            continue
        case, reference = pair
        assert len(case) == size and len(reference) == reference_n
        assert len(set(case)) == size and len(set(reference)) == reference_n and not set(case) & set(reference)
    if design == 'shared_disjoint':
        assert len(set(np.concatenate([pair[0] for pair in pairs]))) == sum(sizes)
        assert all(np.array_equal(pair[1], pairs[0][1]) for pair in pairs)


def test_full_family_bh_and_plus_one_empirical_p():
    p = [.0001, .01, .1, 1, 1]
    np.testing.assert_allclose(core.bh(p), multipletests(p, method='fdr_bh')[1])
    assert core.empirical_p(np.array([0, 1, 1000]), 1000).tolist() == [1 / 1001, 2 / 1001, 1]
    assert core.bh([.01, 1])[0] == .02
    with pytest.raises(ValueError):
        core.bh([.01, np.nan])


@pytest.mark.parametrize('tied,n,m', [(False, 5, 20), (False, 12, 30), (True, 5, 20)])
def test_optimized_u_p_and_g_match_direct_scipy(tied, n, m):
    rng = np.random.default_rng(19)
    x, y = rng.normal(size=(n, 4)), rng.normal(size=(m, 4))
    if tied:
        x, y = np.round(x), np.round(y)
    u, p, g = core.tests_direct(x, y)
    for j in range(4):
        reference = mannwhitneyu(x[:, j], y[:, j], method='auto', alternative='two-sided')
        assert u[j] == reference.statistic
        assert p[j] == pytest.approx(reference.pvalue)
    df = n + m - 2
    correction = gamma(df / 2) / (np.sqrt(df / 2) * gamma((df - 1) / 2))
    sd = np.sqrt(((n - 1) * x.var(0, ddof=1) + (m - 1) * y.var(0, ddof=1)) / df)
    np.testing.assert_allclose(g, correction * (x.mean(0) - y.mean(0)) / sd)
    np.testing.assert_allclose(core.tests_direct(y, x)[2], -g)


def test_missing_data_and_degenerate_variance():
    with pytest.raises(ValueError, match='finite'):
        core.tests_direct(np.array([[np.nan], [1.]]), np.ones((20, 1)))
    _, p, g = core.tests_direct(np.ones((5, 2)), np.ones((20, 2)))
    assert (p == 1).all() and np.isnan(g).all()
    with pytest.raises(ValueError):
        core.effect(np.ones((1, 2)), np.ones((20, 2)))


def test_reference_bank_matches_direct_and_preserves_reference_counts():
    study = synthetic_study()
    bank = ReferenceBank(study)
    selections = {race: np.arange(20) for race in core.RACES}
    p, g = bank.evaluate(selections)
    for race, data in study.races.items():
        for cell_index in data['rows']:
            name = study.cells.disease.iloc[cell_index]
            x = data['x'][data['patients'].disease.eq(name)]
            _, direct_p, direct_g = core.tests_direct(x, data['y'][:20])
            where = slice(cell_index * 4, (cell_index + 1) * 4)
            np.testing.assert_allclose(p[where], direct_p)
            np.testing.assert_allclose(g[where], direct_g)
    selections['Asian'] = np.zeros(20, int)
    with pytest.raises(ValueError, match='replacement'):
        bank.evaluate(selections)


def test_multiclass_permutations_preserve_race_and_block_counts():
    labels = np.array(['A', 'A', 'B', 'C', 'B', 'C'])
    blocks = [np.array([0, 2, 4]), np.array([1, 3, 5])]
    for blocked in [None, blocks]:
        order = core.permuted_order(len(labels), np.random.default_rng(57), blocked)
        assert sorted(labels[order]) == sorted(labels)
        if blocked is not None:
            for block in blocks:
                assert sorted(labels[order[block]]) == sorted(labels[block])
    study = synthetic_study()
    counts = np.zeros(study.total, int)
    permutation_counts(study, 17, 0, False, counts)
    expected = np.zeros(study.total, int)
    for i, (race, data) in enumerate(study.races.items()):
        order = core.permuted_order(len(data['x']), core.rng_for(17, 2, 0, i))
        labels = data['patients'].disease.to_numpy()
        permuted_labels = np.empty_like(labels)
        permuted_labels[order] = labels
        for j, cell_index in enumerate(data['rows']):
            take = permuted_labels == study.cells.disease.iloc[cell_index]
            u, _, _ = core.tests_direct(data['x'][take], data['x'][~take])
            center = take.sum() * (~take).sum() / 2
            expected[cell_index * 4:(cell_index + 1) * 4] = np.abs(u - center) >= np.abs(data['observed_u'][j] - center) - 1e-10
    np.testing.assert_array_equal(counts, expected)


def test_block_center_and_immobile_assignment():
    study = synthetic_study()
    for data in study.races.values():
        expected = np.zeros_like(data['block_center'])
        patients = data['patients']
        bands = pd.cut(patients.conditional_age_years, [3, 5, 7, 10], right=False).astype(str)
        labels = patients.sex_label + '|' + bands
        for label in labels.unique():
            take = labels.eq(label).to_numpy()
            for j, pos in enumerate(data['positions']):
                expected[j] += np.sum(patients.disease.to_numpy()[take] == data['labels'][pos]) * data['ranks'][take].mean(0)
        n = data['ns'][:, None]
        np.testing.assert_allclose(data['block_center'], expected - n * (n + 1) / 2)
    order = core.permuted_order(4, np.random.default_rng(1), [])
    np.testing.assert_array_equal(order, np.arange(4))


def test_seed_reproducibility_checkpoint_resume_and_corruption(tmp_path):
    study = synthetic_study()
    plan = {'seed': 71}
    fresh, resumed = tmp_path / 'fresh', tmp_path / 'resumed'
    fresh.mkdir()
    resumed.mkdir()
    expected, _ = run_reference(study, plan, fresh, 12, False)
    run_reference(study, plan, resumed, 5, False)
    actual, _ = run_reference(study, plan, resumed, 12, True)
    pd.testing.assert_frame_equal(expected, actual)
    expected_permutation = run_permutation(study, plan, fresh, 12, False)
    run_permutation(study, plan, resumed, 5, False)
    actual_permutation = run_permutation(study, plan, resumed, 12, True)
    pd.testing.assert_frame_equal(expected_permutation, actual_permutation)
    expected_placebo = run_placebo(study, plan, fresh, 12, False)
    run_placebo(study, plan, resumed, 5, False)
    actual_placebo = run_placebo(study, plan, resumed, 12, True)
    pd.testing.assert_frame_equal(expected_placebo, actual_placebo)
    state = core.Checkpoints(resumed, 'permutation', core.signature(plan, 'permutation'), True)
    with (state.directory / state.state['chunks'][0]['file']).open('ab') as handle:
        handle.write(b'CORRUPTION')
    with pytest.raises(ValueError, match='corruption'):
        core.Checkpoints(resumed, 'permutation', core.signature(plan, 'permutation'), True)


def test_fabricated_null_and_shifted_positive_control():
    rng = np.random.default_rng(516)
    discoveries, shifted = [], []
    for _ in range(200):
        reference = rng.normal(size=(60, 4))
        cases = rng.normal(size=(20, 4))
        discoveries.append((core.bh(core.tests_direct(cases, reference)[1]) < .05).any())
        shifted.append((core.bh(core.tests_direct(cases + 2, reference)[1]) < .05).sum())
    assert np.mean(discoveries) < .12
    assert np.mean(shifted) > 3.5


def test_missing_baseline_no_fallback_and_patient_independence(tmp_path):
    with pytest.raises(FileNotFoundError, match='BLOCKED: new race-stratified baseline unavailable'):
        core.load_baseline(tmp_path / 'unavailable')
    study = synthetic_study()
    duplicated = pd.concat([study.disease, study.disease.iloc[[0]]])
    with pytest.raises(ValueError, match='Unique'):
        core.baseline.validate_patients(duplicated)
    with pytest.raises(ValueError, match='local_review'):
        core.baseline.private_output(core.ROOT / 'stats_results/unsafe_sensitivity')
    for path in HERE.rglob('*'):
        if path.is_file() and '__pycache__' not in path.parts:
            assert path.suffix in {'.py', '.md'}
    for path in HERE.glob('*.py'):
        assert 'import sklearn' not in path.read_text()
