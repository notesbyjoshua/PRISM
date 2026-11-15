import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.special import gamma
from scipy.stats import mannwhitneyu
from sklearn.impute import SimpleImputer
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import common
import run_five_syndrome_ml as ml
import run_race_stratified as statistics


@pytest.mark.parametrize('source,label,expected', [
    ('GMDB', 'African', 'African'), ('GMDB', 'European', 'European'), ('GMDB', 'Asian', 'Asian'),
    ('FairFace', 'Black', 'African'), ('FairFace', 'White', 'European'),
    ('FairFace', 'East Asian', 'Asian'), ('FairFace', 'Southeast Asian', 'Asian'),
    ('FairFace', 'Indian', 'Asian'), ('FairFace', 'Middle Eastern', 'Asian'),
    ('FairFace', 'West Asian', 'Asian'), ('FairFace', 'Central Asian', 'Asian'),
    ('FairFace', 'Latino_Hispanic', None), ('GMDB', 'Others', None), ('GMDB', 'Unknown', None),
    ('FairFace', 'Other', None), ('GMDB', None, None), ('FairFace', 'unrecognized', None),
])
def test_exact_crosswalk_and_exclusions(source, label, expected):
    assert common.map_group(label, source) == expected


def synthetic_patients(n=5):
    return pd.DataFrame({'patient_id': [f'synthetic-p{i}' for i in range(n)],
                         'image_id': [f'synthetic-image{i}' for i in range(n)],
                         'disease': ['SYNTHETIC syndrome'] * n, 'canonical_race': ['Asian'] * n,
                         'image_dataset_ethnicity_category': ['Asian'] * n})


def test_repeated_patients_or_images_cannot_create_eligibility():
    frame = synthetic_patients(4)
    assert not common.cohort_counts(frame).eligible.any()
    with pytest.raises(ValueError, match='Unique'):
        common.cohort_counts(pd.concat([frame, frame.iloc[[0]]]))
    frame = synthetic_patients()
    frame.loc[4, 'image_id'] = frame.image_id.iloc[0]
    with pytest.raises(ValueError, match='image_id'):
        common.cohort_counts(frame)


def test_counts_include_all_syndromes_and_zero_cells():
    frame = synthetic_patients()
    counts = common.cohort_counts(frame)
    assert len(counts) == 3 and counts.eligible.sum() == 1
    assert counts.n_current_gmdb.sum() == 5
    assert not counts.previous_study_ge10.any()


def test_reference_content_qc_excludes_conflicting_labels_before_deduplication():
    frame = pd.DataFrame({'file': ['SYNTHETIC/b', 'SYNTHETIC/a', 'SYNTHETIC/c', 'SYNTHETIC/d'],
                          'image_sha256': ['synthetic_hash_1', 'synthetic_hash_1', 'synthetic_hash_2', 'synthetic_hash_2'],
                          'race': ['White', 'White', 'Black', 'East Asian']})
    kept, audit = common.screen_reference_content(frame)
    assert kept.file.tolist() == ['SYNTHETIC/a']
    assert audit == {'fairface_conflicting_label_hash_sets': 1, 'fairface_conflicting_label_rows_excluded': 2,
                     'fairface_repeated_content_rows_removed': 1}


def test_quota_capacity_exact_maximality_and_ties():
    composition = {'African': 43, 'Asian': 157, 'European': 536}
    capacity = {'African': 1174, 'Asian': 3227, 'European': 647}
    size, quotas = common.apportion(composition, capacity)
    assert size < 3409 and quotas['European'] == 647 and sum(quotas.values()) == size
    assert common.apportion(composition, capacity) == (size, quotas)
    for candidate in range(size + 1, 3410):
        _, candidate_quotas = common.apportion(composition, {r: 10000 for r in common.RACES}, candidate)
        assert any(candidate_quotas[r] > capacity[r] for r in common.RACES)
    assert common.apportion(dict.fromkeys(common.RACES, 1), dict.fromkeys(common.RACES, 4), 2)[1] == {
        'African': 1, 'Asian': 1, 'European': 0}
    with pytest.raises(ValueError):
        common.apportion(composition, dict.fromkeys(common.RACES, 0))


def test_sampling_is_repeatable_without_replacement_and_matching_is_within_category():
    reference = pd.DataFrame([{'file': f'SYNTHETIC/{race}/{i}', 'canonical_race': race}
                              for race in common.RACES for i in range(30)])
    quotas = {'African': 10, 'Asian': 11, 'European': 12}
    first = common.sample_reference(reference, quotas, 71)
    second = common.sample_reference(reference.sample(frac=1, random_state=2), quotas, 71)
    pd.testing.assert_frame_equal(first, second)
    assert first.file.is_unique and first.canonical_race.value_counts().to_dict() == quotas
    cases, controls = common.matched_frames(synthetic_patients(), first, 'SYNTHETIC syndrome', 'Asian')
    assert len(cases) == 5 and len(controls) == 11 and controls.canonical_race.eq('Asian').all()
    with pytest.raises(ValueError, match='capacity'):
        common.sample_reference(reference, {'African': 31, 'Asian': 1, 'European': 1}, 1)


def test_exact_120_feature_order_and_metadata_exclusion():
    names = [f'synthetic_f{i:03d}' for i in range(120)]
    frame = pd.DataFrame(np.ones((5, 120)), columns=names)
    frame['patient_id'] = 'SYNTHETIC'
    assert common.validate_features(frame, names).columns.tolist() == names
    assert 'patient_id' not in common.validate_features(frame, names)
    with pytest.raises(ValueError, match='120'):
        common.validate_features(frame, names[:-1])
    with pytest.raises(ValueError, match='order'):
        common.validate_features(frame, names[::-1])
    frame.loc[0, names[0]] = np.nan
    with pytest.raises(ValueError, match='Nonfinite'):
        common.validate_features(frame, names)


def test_mannwhitney_and_exact_hedges_sign_and_correction():
    x, y = np.array([4., 5., 6., 8., 9.]), np.arange(10.) / 4
    row = statistics.two_sample(x, y)
    reverse = statistics.two_sample(y, x)
    df = len(x) + len(y) - 2
    correction = gamma(df / 2) / (np.sqrt(df / 2) * gamma((df - 1) / 2))
    sd = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df)
    assert row['hedges_g'] == pytest.approx(correction * (x.mean() - y.mean()) / sd)
    assert row['hedges_g'] > 0 and reverse['hedges_g'] == pytest.approx(-row['hedges_g'])
    assert reverse['p'] == pytest.approx(row['p'])
    assert row['mannwhitney_u'] + reverse['mannwhitney_u'] == len(x) * len(y)
    tied = statistics.two_sample([1, 1, 2, 3, 4], [1, 2, 2, 3, 4, 5])
    assert tied['has_ties'] and tied['mw_method'] == 'asymptotic'
    assert tied['p'] == mannwhitneyu([1, 1, 2, 3, 4], [1, 2, 2, 3, 4, 5], method='asymptotic').pvalue
    untied = statistics.two_sample(np.arange(5.), np.arange(20.) + .1)
    assert untied['mw_method'] == 'exact'
    constant = statistics.two_sample(np.ones(5), np.ones(30))
    assert constant['p'] == 1 and np.isnan(constant['hedges_g'])


def test_planned_hypotheses_and_global_bh_include_unsupported_cells(tmp_path):
    features = [f'synthetic_f{i:03d}' for i in range(120)]
    cases = pd.concat([synthetic_patients(), pd.DataFrame({feature: np.arange(5.) for feature in features})], axis=1)
    reference = pd.DataFrame({feature: np.arange(10.) for feature in features})
    reference['canonical_race'] = 'Asian'
    definitions = pd.DataFrame({'feature': features, 'formula': 'SYNTHETIC', 'units': 'unitless',
                                'anatomical_region': 'SYNTHETIC'})
    result = statistics.primary_tests(cases, reference, common.cohort_counts(cases), features, definitions)
    assert len(result) == 120 and result.p.eq(1).all() and result.q.eq(1).all()
    assert not result.significant.any()
    p = np.array([.001, .02, .03, .9, 1.])
    np.testing.assert_allclose(statistics.bh(p), multipletests(p, method='fdr_bh')[1])
    path = tmp_path / 'synthetic_tests.csv'
    pd.DataFrame({'p': p, 'q': statistics.bh(p)}).to_csv(path, index=False)
    stored = pd.read_csv(path)
    np.testing.assert_allclose(statistics.bh(stored.p), stored.q)
    with pytest.raises(ValueError):
        statistics.bh([.01, np.nan])


def test_analysis_a_has_120_family_and_constant_null():
    features = [f'synthetic_f{i:03d}' for i in range(120)]
    frame = pd.DataFrame(np.ones((60, 120)), columns=features)
    frame['canonical_race'] = np.repeat(common.RACES, 20)
    definitions = pd.DataFrame({'feature': features, 'formula': 'SYNTHETIC', 'units': 'unitless',
                                'anatomical_region': 'SYNTHETIC'})
    result = statistics.analysis_a(frame, features, definitions)
    assert len(result) == 120 and result.q.eq(1).all() and not result.significant.any()


def test_ml_selection_uses_five_different_syndromes_in_one_category_only():
    counts = pd.DataFrame([{'disease': f'SYNTHETIC syndrome {i}', 'canonical_race': race,
                            'n_current_gmdb': (30 - i if race == 'European' else 50 if i == 0 else 9),
                            'eligible': True, 'significant_results': 999 if race == 'Asian' else 0}
                           for race in common.RACES for i in range(6)])
    global_top, selected = common.select_ml(counts)
    assert global_top.canonical_race.nunique() > 1
    assert len(selected) == selected.disease.nunique() == 5 and selected.canonical_race.eq('European').all()
    modified = counts.copy()
    modified['significant_results'] = -100
    pd.testing.assert_frame_equal(selected.drop(columns='significant_results'),
                                  common.select_ml(modified)[1].drop(columns='significant_results'))
    counts['n_current_gmdb'] = 9
    assert common.select_ml(counts)[1].empty


def test_patient_fold_independence_and_training_only_preprocessing(monkeypatch):
    rng = np.random.default_rng(41)
    X = pd.DataFrame(rng.normal(size=(50, 120)), columns=[f'synthetic_f{i:03d}' for i in range(120)])
    y = np.repeat(np.arange(5), 10)
    patients = np.array([f'SYNTHETIC patient {i}' for i in range(50)])
    fit_indices = []

    class RecordingImputer(SimpleImputer):
        def fit(self, X, y=None):
            fit_indices.append(set(X.index))
            return super().fit(X, y)

    monkeypatch.setattr(ml, 'SimpleImputer', RecordingImputer)
    original_spec = ml.model_spec

    def small_grid(name, seed):
        model, _ = original_spec(name, seed)
        return model, {'model__C': [.1]}

    monkeypatch.setattr(ml, 'model_spec', small_grid)
    seen = []

    def observer(estimator, train, test):
        assert len(fit_indices) == 4
        assert all(indices <= set(train) and not indices & set(test) for indices in fit_indices)
        assert all(len(indices) < len(train) for indices in fit_indices[:3])
        assert fit_indices[-1] == set(train)
        np.testing.assert_allclose(estimator.named_steps['imputer'].statistics_, np.median(X.iloc[train], axis=0))
        np.testing.assert_allclose(estimator.named_steps['scaler'].mean_, X.iloc[train].mean())
        seen.extend(test)
        fit_indices.clear()

    result = ml.nested_cv(X, y, patients, 'logistic', observer=observer)
    assert sorted(seen) == list(range(50))
    assert result[4].patient_id.is_unique and len(result[4]) == 50
    assert result[3].filter(regex='test_class_.*_n').ge(1).all().all()
    with pytest.raises(ValueError, match='one row'):
        ml.patient_splits(y, ['same'] * 50, 5, 1)


def test_private_output_guard_and_tracked_artifact_privacy():
    with pytest.raises(ValueError, match='local_review'):
        common.private_output(common.ROOT / 'stats_results/unsafe_redesign')
    tracked = subprocess.check_output(['git', 'ls-files', '--', str(HERE)], cwd=common.ROOT, text=True).splitlines()
    for name in tracked:
        assert Path(name).suffix in {'.py', '.md', '.txt'}
        assert not Path(name).name.startswith('PRIVATE_')
    for path in HERE.rglob('*'):
        if path.is_file() and not any(part in {'__pycache__', '.pytest_cache'} for part in path.parts):
            assert path.suffix in {'.py', '.md', '.txt'}
    assert subprocess.run(['git', 'check-ignore', '-q', str(common.PRIVATE / 'PRIVATE_oof.csv')],
                          cwd=common.ROOT).returncode == 0
