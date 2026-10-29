import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm, rankdata

import extended_analysis as e
import followup as f


def synthetic(n=30):
    rng = np.random.default_rng(91)
    target = np.arange(n*2) < n
    age = rng.uniform(3, 10, size=2*n)
    female = np.tile([0, 1], n)
    return pd.DataFrame({'patient_id': [str(i) for i in range(2*n)],
                         'disease': np.where(target, 'A', 'B'), 'sex_label': np.where(female, 'Female', 'Male'),
                         'conditional_age_years': age, 'x': .3*target+.1*age+.2*female+rng.normal(0, 1, 2*n)})


def manual_hc3(frame, syndrome, feature):
    part = frame.assign(female=frame.sex_label.map({'Male': 0., 'Female': 1.})).dropna(
        subset=[feature, 'conditional_age_years', 'female'])
    x = np.column_stack([np.ones(len(part)), part.disease.eq(syndrome), part.conditional_age_years-6, part.female])
    y = part[feature].to_numpy(float)
    inv = np.linalg.inv(x.T@x)
    beta = inv@x.T@y
    residual = y-x@beta
    h = np.sum((x@inv)*x, axis=1)
    meat = x.T@((residual/(1-h))[:, None]**2*x)
    se = np.sqrt((inv@meat@inv)[1, 1])
    return beta[1], se, 2*norm.sf(abs(beta[1]/se))


def manual_bh(values):
    p = np.asarray(values)
    order = np.argsort(p)
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1])
    return result


def manual_rank_p(x, y):
    combined = np.r_[x, y]
    nx, ny = len(x), len(y)
    u = rankdata(combined)[:nx].sum()-nx*(nx+1)/2
    ties = np.unique(combined, return_counts=True)[1]
    variance = nx*ny/12*(nx+ny+1-np.sum(ties**3-ties)/((nx+ny)*(nx+ny-1)))
    p = min(1., 2*norm.sf((abs(u-nx*ny/2)-.5)/np.sqrt(variance))) if variance > 0 else 1.
    return u, p


def direct_g(x, y):
    nx, ny = len(x), len(y)
    pooled = ((nx-1)*np.var(x, ddof=1)+(ny-1)*np.var(y, ddof=1))/(nx+ny-2)
    return (1-3/(4*(nx+ny)-9))*(np.mean(x)-np.mean(y))/np.sqrt(pooled)


def explicit_bootstrap(x, y, count, seed):
    rng = np.random.default_rng(seed)
    a = x[rng.integers(len(x), size=(count, len(x))) ]
    b = y[rng.integers(len(y), size=(count, len(y))) ]
    variance = ((len(x)-1)*a.var(axis=1, ddof=1)+(len(y)-1)*b.var(axis=1, ddof=1))/(len(x)+len(y)-2)
    draws = (1-3/(4*(len(x)+len(y))-9))*(a.mean(axis=1)-b.mean(axis=1))/np.sqrt(variance)
    return np.quantile(draws[np.isfinite(draws)], [.025, .975])


def test_pair_family_is_unique_and_includes_all_six_anchor_pairs():
    targets = [*e.ANCHORS, *[f'other_{i}' for i in range(45)]]
    pairs = e.planned_pairs(targets)
    assert len(pairs) == 186
    assert len({frozenset(p) for p in pairs}) == 186
    assert sum(a in e.ANCHORS and b in e.ANCHORS for a, b in pairs) == 6
    assert all(a != b and (a in e.ANCHORS or b in e.ANCHORS) for a, b in pairs)


def test_unknown_anchor_is_not_silently_dropped():
    with pytest.raises(ValueError, match='anchor'):
        e.planned_pairs(['A', 'B'])


def test_complete_adjusted_model_matches_independent_hc3():
    frame = synthetic()
    row = e.adjusted_fit(frame, 'A', 'x')
    beta, se, p = manual_hc3(frame, 'A', 'x')
    assert row['status'] == 'ok'
    np.testing.assert_allclose([row['target_beta'], row['hc3_se'], row['p_value']], [beta, se, p], atol=1e-12)
    np.testing.assert_allclose([row['ci_95_lower'], row['ci_95_upper']], beta+np.array([-1, 1])*norm.ppf(.975)*se)


def test_missing_covariates_counts_and_exclusion_support():
    frame = synthetic()
    frame.loc[:10, 'sex_label'] = None
    row = e.adjusted_fit(frame, 'A', 'x')
    assert row['target_missing_covariates'] == 11 and row['n_target'] == 19
    assert row['status'] == 'unsupported' and 'exclusion' in row['reason']
    assert np.isnan(row['target_beta'])


def test_nonoverlap_is_not_fitted():
    frame = synthetic()
    frame.loc[frame.disease.eq('A'), 'conditional_age_years'] = np.linspace(3, 5, 30)
    frame.loc[frame.disease.eq('B'), 'conditional_age_years'] = np.linspace(6, 9, 30)
    row = e.adjusted_fit(frame, 'A', 'x')
    assert row['status'] == 'unsupported' and 'age_overlap' in row['reason']


def test_rank_deficiency_is_flagged():
    frame = synthetic()
    frame['sex_label'] = np.where(frame.disease.eq('A'), 'Female', 'Male')
    row = e.adjusted_fit(frame, 'A', 'x')
    assert row['status'] == 'unsupported' and 'rank' in row['reason']


def test_constant_outcome_is_not_inferred():
    frame = synthetic()
    frame['x'] = 1.
    row = e.adjusted_fit(frame, 'A', 'x')
    assert row['status'] == 'unsupported' and 'residual' in row['reason']


def test_correction_keeps_unsupported_tests_in_family():
    frame = pd.DataFrame({'status': ['ok', 'unsupported', 'ok'], 'p_value': [.01, np.nan, .04]})
    result = e.corrected(frame, 'fixture', 3)
    np.testing.assert_allclose(result.correction_p, [.01, 1, .04])
    np.testing.assert_allclose(result.q_global_bh.dropna(), [.03, .06])
    assert not result.significant.iloc[1]
    with pytest.raises(ValueError):
        e.corrected(frame, 'fixture', 2)


def test_duplicate_recorded_ids_rejected():
    frame = synthetic()
    frame.loc[1, 'patient_id'] = frame.loc[0, 'patient_id']
    with pytest.raises(ValueError, match='IDs'):
        e.require_unique(frame)


def test_pairwise_overlap_rejected():
    frame = synthetic()
    with pytest.raises(ValueError, match='overlap'):
        e.two_groups(frame.iloc[:30], frame.iloc[20:], ['x'], 100, 1)


def test_two_groups_rank_g_and_bootstrap_are_correct():
    frame = synthetic()
    a, b = frame.iloc[:30], frame.iloc[30:]
    row = e.two_groups(a, b, ['x'], 100, 9)[0]
    x, y = a.x.to_numpy(), b.x.to_numpy()
    u, p = manual_rank_p(x, y)
    assert row['status'] == 'ok'
    np.testing.assert_allclose([row['hedges_g'], row['mannwhitney_u'], row['p_value']], [direct_g(x, y), u, p])
    np.testing.assert_allclose([row['ci_95_lower'], row['ci_95_upper']], explicit_bootstrap(x, y, 100, 9))


def test_sparse_stratum_exports_status_and_no_fabricated_ci():
    frame = synthetic()
    row = e.two_groups(frame.iloc[:9], frame.iloc[30:], ['x'], 100, 1, min_b=20)[0]
    assert row['status'] == 'unsupported' and row['bootstrap_finite'] == 0
    assert np.isnan(row['hedges_g']) and np.isnan(row['ci_95_lower'])
    empty = e.two_groups(frame.iloc[:0], frame.iloc[30:], ['x'], 100, 1)[0]
    assert empty['status'] == 'unsupported'


def test_shared_moments_reproduce_direct_contrasts():
    rng = np.random.default_rng(10)
    values, ref = rng.normal(size=(60, 3)), rng.normal(size=(50, 3))
    fair, within = e.shared_effects(values, [0, 20, 40], np.array([0, 1]), ref, np.array([0, 1, 2]))
    for i in range(2):
        target = values[20*i:20*(i+1)]
        other = np.delete(values, np.s_[20*i:20*(i+1)], axis=0)
        np.testing.assert_allclose(fair[i], f.effect(target, ref), atol=1e-12)
        np.testing.assert_allclose(within[i], f.effect(target, other), atol=1e-12)
    _, restricted = e.shared_effects(values, [0, 20, 40], np.array([0, 1]), ref, np.array([0, 1]))
    np.testing.assert_allclose(restricted[0], -restricted[1], atol=1e-12)


def test_attenuation_has_explicit_direction():
    x, y = np.array([1., -2., .5]), np.array([.5, -1., .2])
    values = dict(zip(e.METRICS, e.agreement(x, y)))
    assert values['direction_agreement'] == 1
    assert values['fraction_absolute_attenuated'] == 1
    assert values['median_absolute_attenuation'] == .5


def test_failed_input_is_logged_without_masking_exception(tmp_path, monkeypatch):
    source = tmp_path/'review'
    source.mkdir()
    monkeypatch.setattr(f, 'OUT', source)
    monkeypatch.setattr(f, 'facekit_provenance', lambda: {})
    monkeypatch.setattr(e, 'version_record', lambda: {})
    monkeypatch.setattr(e.subprocess, 'check_output', lambda *a, **kw: 'fixture')
    monkeypatch.setattr('sys.argv', ['extended_analysis.py', '--output-dir', str(source/'local_only/failure')])

    def missing_load():
        f.INPUTS.add(f.ROOT/'absent.csv')
        raise FileNotFoundError('missing fixture input')

    monkeypatch.setattr(f, 'load', missing_load)
    with pytest.raises(FileNotFoundError, match='missing fixture input'):
        e.main()
    out = source/'local_only/failure'
    status = json.loads((out/'execution_status.json').read_text())
    provenance = json.loads((out/'provenance.json').read_text())
    assert status['status'] == 'failed' and status['exit_code'] == 1
    assert provenance['missing_inputs'] == ['absent.csv']


def test_output_refuses_overwrite_and_outside_local_only(tmp_path, monkeypatch):
    monkeypatch.setattr(f, 'OUT', tmp_path)
    existing = tmp_path/'local_only/existing'
    existing.mkdir(parents=True)
    (existing/'marker').write_text('unchanged')
    monkeypatch.setattr('sys.argv', ['extended_analysis.py', '--output-dir', str(existing)])
    with pytest.raises(ValueError, match='overwrite'):
        e.main()
    assert list(existing.iterdir()) == [existing/'marker']
    outside = tmp_path/'outside'
    monkeypatch.setattr('sys.argv', ['extended_analysis.py', '--output-dir', str(outside)])
    with pytest.raises(ValueError, match='local_only'):
        e.main()
    assert not outside.exists()
