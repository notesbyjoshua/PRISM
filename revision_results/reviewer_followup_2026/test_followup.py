import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import followup as f


@pytest.fixture(scope='module')
def cohort():
    return f.load()


def test_comparator_excludes_target_and_patients():
    frame = pd.DataFrame({'disease': ['A', 'A', 'B', 'C'], 'patient_id': ['1', '2', '3', '4']})
    target, other = f.comparator_frames(frame, 'A', ['A', 'B'], 'all')
    assert target.patient_id.tolist() == ['1', '2']
    assert other.patient_id.tolist() == ['3', '4']
    _, restricted = f.comparator_frames(frame, 'A', ['A', 'B'], 'eligible49')
    assert restricted.patient_id.tolist() == ['3']


def test_unique_selected_patients(cohort):
    raw, _, primary, alternate, _, _, _, _ = cohort
    assert primary.patient_id.is_unique and primary.patient_id.notna().all()
    assert set(primary.patient_id) == set(alternate.patient_id)
    assert f.clean_id(raw.patient_id).isna().sum() == 0
    assert len(primary) == 1850


def test_exact_historical_features(cohort):
    features = cohort[-1]
    assert len(features) == len(set(features)) == 125
    assert cohort[2].loc[:, features[0]:features[-1]].columns.tolist() == features
    assert set(features) == set(cohort[5].feature)


def test_hedges_sign_and_hand_calculation():
    x, y = np.array([5., 6., 7.]), np.array([1., 2., 3., 4.])
    expected = (1-3/(4*7-9))*(6-2.5)/np.sqrt((2*1+3*(5/3))/5)
    assert np.isclose(f.effect(x, y), expected)
    assert np.isclose(f.effect(y, x), -expected)


def test_bh_manual_and_full_family():
    p = np.array([.04, .001, .2, .015])
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order]*4/np.arange(1, 5))[::-1])[::-1]
    expected = np.empty(4)
    expected[order] = np.minimum(q, 1)
    assert np.allclose(f.bh(p), expected)
    result = pd.read_csv(f.OUT/'within_gmdb_results.csv')
    assert len(result) == 6125 and not result.duplicated(['syndrome', 'feature']).any()
    assert np.allclose(f.bh(result.p_value), result.q_global_bh)


def test_bootstrap_is_patient_row_resampling():
    values = np.column_stack([np.arange(1., 8.), np.arange(1., 8.)*3])
    rng1, rng2 = np.random.default_rng(123), np.random.default_rng(123)
    count, mean, var = f.bootstrap_moments(values, 30, rng1)
    indexes = rng2.integers(len(values), size=(30, len(values)))
    expected = values[indexes]
    assert np.array_equal(count, np.full((30, 2), 7))
    assert np.allclose(mean, expected.mean(axis=1))
    assert np.allclose(var, expected.var(axis=1, ddof=1))
    assert np.allclose(mean[:, 1], 3*mean[:, 0])


def test_bootstrap_independent_scalar_and_seed():
    x = np.arange(1., 12.)[:, None]
    y = np.arange(2., 19.)[:, None]
    actual = f.bootstrap_effect(x, y, 80, [123, 9])
    rng = np.random.default_rng([123, 9])
    xi = rng.integers(len(x), size=(80, len(x)))
    yi = rng.integers(len(y), size=(80, len(y)))
    expected = np.array([f.effect(x[a], y[b]) for a, b in zip(xi, yi)])
    assert np.allclose(actual, expected, rtol=1e-10, atol=1e-10)
    assert np.array_equal(actual, f.bootstrap_effect(x, y, 80, [123, 9]))
    assert not np.array_equal(actual, f.bootstrap_effect(x, y, 80, [123, 10]))


def test_bootstrap_missing_values_are_featurewise():
    x = np.array([[1., np.nan], [2., 4.], [3., 6.], [4., 8.]])
    count, mean, var = f.bootstrap_moments(x, 15, np.random.default_rng(7))
    assert (count[:, 0] == 4).all()
    assert (count[:, 1] <= 4).all()
    assert np.isfinite(mean[:, 0]).all()


def test_independent_subset_and_frozen_baseline(cohort):
    _, _, primary, _, ref, frozen, targets, features = cohort
    result = pd.read_csv(f.OUT/'within_gmdb_results.csv')
    validations = []
    selected = [(targets[0], features[0]), (targets[13], features[29]), (targets[-1], features[-1])]+f.HIGHLIGHTS
    for syndrome, feature in selected:
        target, other = f.comparator_frames(primary, syndrome, targets, 'all')
        x, y = target[feature].to_numpy(float), other[feature].to_numpy(float)
        n, m = len(x), len(y)
        pooled = ((n-1)*x.var(ddof=1)+(m-1)*y.var(ddof=1))/(n+m-2)
        g = (1-3/(4*(n+m)-9))*(x.mean()-y.mean())/np.sqrt(pooled)
        z = np.r_[x, y]
        ranks = rankdata(z)
        u = ranks[:n].sum()-n*(n+1)/2
        _, ties = np.unique(z, return_counts=True)
        tie_sum = (ties**3-ties).sum()
        sd = np.sqrt(n*m/12*((n+m+1)-tie_sum/((n+m)*(n+m-1))))
        p = min(1., 2*norm.sf((abs(u-n*m/2)-.5)/sd))
        row = result.loc[result.syndrome.eq(syndrome) & result.feature.eq(feature)].iloc[0]
        assert np.isclose(g, row.hedges_g, rtol=1e-10)
        assert np.isclose(p, row.p_value, rtol=1e-8, atol=1e-14)
        old = frozen.loc[frozen.disease.eq(syndrome) & frozen.feature.eq(feature)].iloc[0]
        assert np.isclose(f.effect(x, ref[feature].to_numpy(float)), old.hedges_g, rtol=1e-10)
        validations.append({'syndrome': syndrome, 'feature': feature, 'independent_g': g,
                            'independent_mannwhitney_p': p, 'matches_new_and_frozen': True})
    f.write(pd.DataFrame(validations), 'independent_subset_validation.csv')


def test_historical_dictionary_rebuilt_from_exact_git_source():
    path = f.ROOT/'publication/scripts/generate_feature_definitions.py'
    spec = importlib.util.spec_from_file_location('historical_definitions', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source, tree, namespace, records = module.load_historical_source()
    schema = pd.read_csv(f.FROZEN/'exact_feature_names.csv').feature.tolist()
    dependencies = module.landmark_dependencies(namespace, schema)
    rows, _ = module.formula_rows(source, tree, namespace, records, dependencies)
    current = pd.DataFrame(rows).fillna('').astype(str)
    stored = pd.read_csv(f.OUT/'feature_dictionary.csv').fillna('').astype(str)
    pd.testing.assert_frame_equal(current[stored.columns], stored, check_dtype=False)
    (f.OUT/'historical_dictionary_validation.json').write_text(__import__('json').dumps(
        {'features': 125, 'all_definition_columns_match': True, 'source_commit': module.COMMIT,
         'historical_source_sha256': __import__('hashlib').sha256(source.encode()).hexdigest()}, indent=2))


def test_independent_hc3_interactions(cohort):
    _, _, primary, _, ref, _, _, _ = cohort
    validations = []
    for layer, directory in f.LAYERS.items():
        saved = pd.read_csv(f.FROZEN/f'age_3_9/ethnicity_analysis/{directory}/analysis_C_all_interactions.csv.gz')
        for row in saved.iloc[[0, len(saved)//2, -1]].itertuples():
            disease = primary.loc[primary.internal_syndrome_name.astype('string').str.strip().eq(row.syndrome_name)]
            vals, status, group = [], [], []
            for is_disease, data in [(0, ref), (1, disease)]:
                for is_group2, name in [(0, row.group_1), (1, row.group_2)]:
                    y = data.loc[data[f'{layer}_group'].eq(name), row.feature].dropna().to_numpy(float)
                    vals.extend(y)
                    status.extend([is_disease]*len(y))
                    group.extend([is_group2]*len(y))
            d, g, y = np.array(status), np.array(group), np.array(vals)
            design = np.column_stack([np.ones(len(d)), d, g, d*g])
            inv = np.linalg.inv(design.T@design)
            beta = inv@design.T@y
            residual = y-design@beta
            leverage = np.einsum('ij,jk,ik->i', design, inv, design)
            meat = design.T@(((residual/(1-leverage))**2)[:, None]*design)
            se = np.sqrt(np.diag(inv@meat@inv))[-1]
            p = 2*norm.sf(abs(beta[-1]/se))
            assert np.isclose(beta[-1], row.interaction_beta_group2_minus_group1, rtol=1e-7, atol=1e-10)
            assert np.isclose(se, row.interaction_standard_error, rtol=1e-7)
            assert np.isclose(p, row.interaction_p_value, rtol=1e-6, atol=1e-12)
            validations.append({'layer': layer, 'syndrome': row.syndrome_name, 'feature': row.feature,
                                'group_1': row.group_1, 'group_2': row.group_2, 'n': len(y), 'passed': True})
    f.write(pd.DataFrame(validations), 'independent_hc3_validation.csv')


def test_no_patient_identifiers_in_aggregate_csvs():
    prohibited = {'patient_id', 'image_id', 'patient_name', 'reference_key', 'facekit_relative_path', 'hgvs', 'hpo_terms'}
    for path in f.OUT.glob('*.csv'):
        assert not prohibited & set(pd.read_csv(path, nrows=0).columns), path.name


def test_all_frozen_primary_effects_pvalues_and_bh(cohort):
    _, _, primary, _, ref, frozen, targets, features = cohort
    rows = []
    for syndrome in targets:
        x = primary.loc[primary.disease.eq(syndrome), features].to_numpy(float)
        y = ref[features].to_numpy(float)
        test = f.mannwhitneyu(x, y, axis=0, method='asymptotic', alternative='two-sided')
        for feature, g, p in zip(features, f.effect(x, y), test.pvalue):
            rows.append({'disease': syndrome, 'feature': feature, 'g': g, 'p': p})
    actual = pd.DataFrame(rows)
    actual['q'] = f.bh(actual.p)
    joined = actual.merge(frozen, on=['disease', 'feature'], validate='one_to_one')
    assert len(joined) == 6125
    assert np.allclose(joined.g, joined.hedges_g, rtol=1e-10, atol=1e-12)
    assert np.allclose(joined.p, joined.mannwhitney_p, rtol=1e-9, atol=1e-14)
    assert np.allclose(joined.q, joined.mannwhitney_q_fdr_global, rtol=1e-9, atol=1e-14)
    assert joined.q.lt(.05).sum() == 811
    assert (joined.q.lt(.05) & joined.g.abs().ge(.8)).sum() == 379
    (f.OUT/'baseline_numerical_validation.json').write_text(__import__('json').dumps(
        {'primary_test_count': 6125, 'all_g_match': True, 'all_p_match': True, 'all_q_match': True,
         'significant': 811, 'significant_abs_g_ge_08': 379}, indent=2))


def test_reference_omnibus_families(cohort):
    from scipy.stats import kruskal
    ref, features = cohort[4], cohort[-1]
    for layer, directory in f.LAYERS.items():
        groups = [data[features].to_numpy(float) for _, data in ref.groupby(f'{layer}_group') if len(data) >= 20]
        actual = kruskal(*groups, axis=0)
        saved = pd.read_csv(f.FROZEN/f'age_3_9/ethnicity_analysis/{directory}/analysis_A_healthy_kruskal_all_groups.csv.gz').set_index('feature').loc[features]
        assert np.allclose(actual.statistic, saved.kruskal_H, rtol=1e-9)
        assert np.allclose(actual.pvalue, saved.p_value, rtol=1e-9, atol=1e-14)
        assert np.allclose(f.bh(actual.pvalue), saved.fdr_q_value, rtol=1e-9, atol=1e-14)
