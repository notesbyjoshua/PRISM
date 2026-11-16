import numpy as np
import pandas as pd
import pytest

from stats_testing.robustness_20261110 import common as c


def test_canonical_join_refuses_duplicate_or_incomplete_grid():
    left = pd.DataFrame({'canonical_syndrome_id': ['1', '1'], 'feature': ['a', 'b']})
    right = left.assign(effect=[.5, -.5])
    merged = c.keyed_merge(left, right, 'other_', ['effect'])
    assert merged.other_effect.tolist() == [.5, -.5]
    with pytest.raises(c.ResearchContractError):
        c.keyed_merge(left, right.iloc[[0]], 'other_', ['effect'])
    with pytest.raises(c.ResearchContractError):
        c.keyed_merge(left, pd.concat([right, right.iloc[[0]]]), 'other_', ['effect'])
    absent = c.keyed_merge(left, right.iloc[[0]], 'optional_', ['effect'], False)
    assert np.isnan(absent.optional_effect.iloc[1])


def test_missing_support_keeps_missing_display_q_and_planned_denominator():
    frame = pd.DataFrame({'p': [.001, np.nan, .04]})
    result = c.corrected(frame, 'p', [True, False, True], 'fixture_planned_three')
    assert result.planned_hypotheses.eq(3).all()
    assert result.correction_p.iloc[1] == 1
    assert np.isnan(result.exploratory_q.iloc[1])
    assert result.exploratory_q.iloc[0] == pytest.approx(.003)
    with pytest.raises(c.ResearchContractError):
        c.bh([.01, np.nan])


def test_no_individual_rows_or_nested_values_exported(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'ROOT', tmp_path)
    destination = tmp_path / 'local_review/test/aggregate.csv'
    for frame in [pd.DataFrame({'patient_id': ['not-a-real-patient']}), pd.DataFrame({'values': [[1, 2, 3]]})]:
        with pytest.raises(c.ResearchContractError):
            c.write_aggregate(frame, destination)
    with pytest.raises(c.ResearchContractError):
        c.write_aggregate(pd.DataFrame({'n': [20]}), tmp_path / 'publication/not-approved.csv')
    c.write_aggregate(pd.DataFrame({'n': [20], 'effect': [.3]}), destination)
    assert pd.read_csv(destination).n.tolist() == [20]


def test_effect_sign_small_sample_and_zero_variance():
    x, y = np.arange(10.)[:, None] + 1, np.arange(10.)[:, None]
    expected = (1 - 3 / 71) / np.std(np.arange(10.), ddof=1)
    assert c.effect(x, y)[0] == pytest.approx(expected)
    assert np.isnan(c.effect(np.ones((10, 1)), np.zeros((10, 1)))[0])
    assert np.isnan(c.effect(x[:1], y)[0])


def test_current_age_rules_and_canonical_taxonomy():
    frame = pd.DataFrame({'age_year': [3, 9, 9, 3.5, None], 'age_month': [0, 11, 12, 0, 0]})
    age = c.age_years(frame)
    assert age.iloc[0] == 3 and age.iloc[1] == pytest.approx(9 + 11 / 12)
    assert age.iloc[2:].isna().all()
    with pytest.raises(c.ResearchContractError):
        c.canonicalize(pd.DataFrame({'syndrome': ['unknown display alias']}))


def test_negative_control_pools_are_disjoint_reproducible_and_not_pitch_ordered():
    from stats_testing.robustness_20261110.source_bias import control_indices
    strata, pitch, sizes = np.repeat('fixture', 200), np.arange(200.), [10] * 8
    groups, reference, imbalance = control_indices(strata, pitch, sizes, np.random.default_rng(3), 'pitch_imbalanced')
    again, again_ref, _ = control_indices(strata, pitch, sizes, np.random.default_rng(3), 'pitch_imbalanced')
    assert all(np.array_equal(a, b) for a, b in zip(groups, again)) and np.array_equal(reference, again_ref)
    assert not np.isin(np.concatenate(groups), reference).any()
    assert len(np.unique(np.concatenate(groups))) == 80
    assert not np.all(np.diff([pitch[group].mean() for group in groups]) > 0)
    assert pitch[np.concatenate(groups)].min() >= 100 and pitch[reference].max() < 100
    assert imbalance == 0


def test_patient_bootstrap_preserves_rows_and_rejects_repeat_subjects():
    from stats_testing.robustness_20261110.stability import patient_samples, fixed_reference_effect
    x = np.column_stack([np.arange(10.), 2 * np.arange(10.)])
    samples = patient_samples(x, np.arange(10), 20, np.random.default_rng(9))
    again = patient_samples(x, np.arange(10), 20, np.random.default_rng(9))
    np.testing.assert_array_equal(samples, again)
    np.testing.assert_array_equal(samples[:, :, 1], 2 * samples[:, :, 0])
    with pytest.raises(c.ResearchContractError):
        patient_samples(x, np.repeat(1, 10), 20, np.random.default_rng(9))
    reference = x - 1
    expected = np.array([c.effect(sample, reference) for sample in samples])
    np.testing.assert_allclose(fixed_reference_effect(samples, reference), expected)


def test_normalization_direction_and_simulated_components():
    from stats_testing.robustness_20261110.geometry import alternative_ratios, numerical_geometry_simulation
    frame = pd.DataFrame({'outer_canthal_distance': [60 / 140], 'philtrum_width': [20 / 140],
                          'face_aspect_ratio': [140 / 200], 'midface_width': [100 / 140], 'mouth_width': [50 / 140]})
    outer = alternative_ratios(frame, 'outer_canthal_distance')
    width = alternative_ratios(frame, 'philtrum_width')
    assert outer['height_normalized'][0] == pytest.approx(60 / 200)
    assert outer['midface_width_normalized'][0] == pytest.approx(.6)
    assert width['mouth_width_normalized'][0] == pytest.approx(.4)
    simulation = numerical_geometry_simulation()
    row = simulation.loc[simulation.simulated_numerator_change.eq(0) & np.isclose(simulation.simulated_denominator_change, .1)].iloc[0]
    assert row.relative_ratio_change == pytest.approx(1 / 1.1 - 1)


def test_dependence_graph_and_cluster_labels():
    from stats_testing.robustness_20261110.dependence import dependency_edges, clusters, adjusted_rand
    definitions = pd.DataFrame({'feature': ['a', 'b', 'mean'],
                                 'formula': ['x/bizygomatic_width', 'y/bizygomatic_width', "0.5*(f['a']+f['b'])"]})
    edges = dependency_edges(definitions)
    assert edges.relation.eq('explicit_algebraic_reference').sum() == 2
    assert edges.relation.eq('shared_normalization_not_an_identity').sum() == 2
    with pytest.raises(c.ResearchContractError):
        dependency_edges(definitions.assign(formula="f['unknown']"))
    labels = clusters(np.array([[1., 1., 0.], [1., 1., 0.], [0., 0., 1.]]))
    assert labels[0] == labels[1] and labels[0] != labels[2]
    assert adjusted_rand([1, 1, 2], [6, 6, 7]) == pytest.approx(1)
    exact = dependency_edges(pd.DataFrame({'feature': ['a'], 'formula': ['x/midface_height']}))
    assert exact.parent.tolist() == ['normalizer:midface_height']
