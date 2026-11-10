import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('synthetic_consistency', ROOT / 'stats_testing/synthetic_consistency.py')
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)


def test_signed_shifts_and_fixed_reference_effects():
    real = np.arange(10, dtype=float)
    synthetic = real + 2
    reference = real - 1
    result = sc.feature_metrics(real, synthetic, reference)
    assert result['mean_shift'] == result['median_shift'] == 2
    assert result['variance_ratio'] == pytest.approx(1)
    assert result['wasserstein'] == pytest.approx(2)
    assert result['synthetic_g'] > result['real_g'] > 0
    assert result['synthetic_reference_sd_contrast'] == pytest.approx(3 / reference.std(ddof=1))
    assert result['real_reference_sd_contrast'] == pytest.approx(1 / reference.std(ddof=1))


def test_variance_scaling_does_not_change_fixed_reference_contrast():
    real = np.arange(10, dtype=float)
    synthetic = real.mean() + .2 * (real - real.mean())
    result = sc.feature_metrics(real, synthetic, real - 1)
    assert result['variance_ratio'] == pytest.approx(.04)
    assert result['synthetic_reference_sd_contrast'] == pytest.approx(result['real_reference_sd_contrast'])
    assert result['synthetic_g'] != pytest.approx(result['real_g'])


@pytest.mark.parametrize('real', [np.ones(10), np.ones(1), np.repeat(np.nan, 10)])
def test_zero_variance_missing_and_small_samples_are_undefined(real):
    row = sc.feature_metrics(real, np.ones(10), np.arange(10.))
    assert not row['real_sd_defined']
    assert np.isnan(row['variance_ratio'])
    assert np.isnan(row['standardized_mean_shift'])
    if len(real) == 1 or not np.isfinite(real).any():
        assert row['small_n']


def test_nonfinite_and_wrong_schema_fail_closed():
    features = [f'f{i}' for i in range(120)]
    frame = pd.DataFrame(np.ones((10, 120)), columns=features)
    sc.validate_schema(frame, features, set())
    for changed in [frame.assign(legacy_feature=1), frame.drop(columns='f0'), frame[features[::-1]],
                    frame.assign(f0=np.nan), frame.rename(columns={'f0': 'renamed'})]:
        with pytest.raises(ValueError):
            sc.validate_schema(changed, features, set())


def real_fixture():
    return pd.DataFrame({'patient_id': ['a', 'b'], 'image_id': ['x', 'y'], 'frontal_ok': [True, True],
                         'conditional_age_years': [3., 9.9]})


def test_repeated_patients_and_synthetic_primary_are_rejected():
    sc.validate_real(real_fixture())
    for frame in [real_fixture().assign(patient_id='a'), real_fixture().assign(patient_id=None),
                  real_fixture().assign(synthetic=True), real_fixture().assign(conditional_age_years=10)]:
        with pytest.raises(ValueError):
            sc.validate_real(frame)


def test_mapping_excludes_unknown_and_ambiguous_and_rejects_multiple():
    mapping = pd.DataFrame([{'pdidb_disease': 'a', 'gmdb_syndrome': 'real-a',
                             'mapping_status': 'reviewed_exact_name', 'review_basis': 'exact'},
                            {'pdidb_disease': 'ambiguous', 'gmdb_syndrome': None,
                             'mapping_status': 'excluded_ambiguous_scope', 'review_basis': 'distinct scopes'}])
    result = sc.reviewed_mapping(mapping, {'real-a'}, {'a', 'ambiguous', 'unknown'})
    assert result.set_index('pdidb_disease').loc['unknown', 'mapping_status'] == 'excluded_unmatched'
    assert result.mapping_status.eq('reviewed_exact_name').sum() == 1
    with pytest.raises(ValueError, match='Multiply mapped'):
        sc.reviewed_mapping(pd.concat([mapping, mapping.iloc[[0]]]), {'real-a'}, {'a'})
    alias = mapping.iloc[[0]].assign(pdidb_disease='b')
    with pytest.raises(ValueError, match='one-to-one'):
        sc.reviewed_mapping(pd.concat([mapping, alias]), {'real-a'}, {'a', 'b'})


def test_constant_profiles_have_undefined_correlations():
    assert all(np.isnan(sc.correlations(np.ones(120), np.arange(120))))
    assert all(np.isnan(sc.correlations([1, 2], [3, 4])))


def test_aggregate_export_rejects_identifiers_and_nested_values(tmp_path):
    for frame in [pd.DataFrame({'image_id': ['private']}), pd.DataFrame({'values': [[1, 2, 3]]})]:
        with pytest.raises(ValueError):
            sc.save_aggregate(frame, tmp_path / 'bad.csv')
    frame = pd.DataFrame({'syndrome': ['aggregate'], 'mean': [1. / 3], 'n': [20]})
    sc.save_aggregate(frame, tmp_path / 'first.csv')
    sc.save_aggregate(frame, tmp_path / 'second.csv')
    assert (tmp_path / 'first.csv').read_bytes() == (tmp_path / 'second.csv').read_bytes()


def test_low_n_model_tags_omitted_and_deterministic():
    features = [f'f{i}' for i in range(120)]
    values = np.arange(1200, dtype=float).reshape(10, 120) / 1000
    real = pd.DataFrame(values, columns=features).assign(disease='real')
    reference = pd.DataFrame(values - 1, columns=features)
    synthetic = pd.DataFrame(np.concatenate([values + .2, values[:2]]), columns=features)
    synthetic['disease'], synthetic['model'] = 'synth', ['v1'] * 10 + ['v2'] * 2
    definitions = pd.DataFrame({'feature': features, 'units': 'unitless', 'formula': 'fixture', 'anatomical_region': 'fixture'})
    mapping = pd.DataFrame({'pdidb_disease': ['synth'], 'gmdb_syndrome': ['real'], 'mapping_status': ['reviewed_exact_name']})
    a, b = sc.analyze(real, reference, synthetic, definitions, mapping)
    again_a, again_b = sc.analyze(real, reference, synthetic, definitions, mapping)
    pd.testing.assert_frame_equal(a, again_a)
    pd.testing.assert_frame_equal(b, again_b)
    assert set(a.model_tag) == {'pooled', 'v1'}
    assert set(a.reference_n) == {10}
    assert b.sensitivity.eq('winsorized_01_99').sum() == 1
