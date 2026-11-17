import importlib.util
import json
from pathlib import Path

import pytest

CODE=Path(__file__).resolve().parents[1]/'generate_prism_sensitivity_figure.py'
spec=importlib.util.spec_from_file_location('sensitivity_figure',CODE)
v=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


@pytest.fixture
def summary():
    return json.loads(v.SOURCE.read_text())


def test_authoritative_summary_and_derived_counts(summary):
    assert v.sha(v.SOURCE)==v.EXPECTED_SOURCE_SHA256
    x=v.values(summary)
    assert (x['zero'],x['nonzero'],x['B'])==(987,13,1000)
    assert (x['unblocked'],x['overlap_unblocked'],x['other_unblocked'])==(103,80,23)
    assert (x['blocked'],x['overlap_blocked'],x['other_blocked'])==(92,58,34)
    assert (x['varied_stable'],x['varied_total'])==(269,387)
    assert x['varied_percent']==pytest.approx(69.5090439276)
    assert x['overall_stable']==849 and x['certain']==0 and x['cross']==55


@pytest.mark.parametrize('change',[
    lambda s:s.update(status='INCOMPLETE'),
    lambda s:s.update(baseline_id='OLD_POOLED_STUDY'),
    lambda s:s['placebo'].update(supported_hypotheses=9600),
    lambda s:s['placebo'].update(fraction_at_least_one=.0135),
    lambda s:s['reference_stability_by_race'][2].update(reference_images_vary=True),
    lambda s:s['reference_stability_by_race'][0].update(descriptive_stable=77),
    lambda s:s['gmdb_permutation']['within_race'].update(overlap_main_same_direction=104),
    lambda s:s['gmdb_permutation']['within_race'].update(MC_upper_q_below_005=1),
    lambda s:s.update(simultaneous_MC_certain_cross_design_candidates=1),
])
def test_inconsistent_source_fails(summary,change):
    change(summary)
    with pytest.raises(ValueError):
        v.validate(summary)


def test_missing_source_field_never_gets_a_default(summary):
    del summary['reference_resampling']['B']
    with pytest.raises(KeyError):
        v.values(summary)


@pytest.fixture
def comparison_source(summary):
    plan = {'baseline_id': summary['baseline_id'], 'planned_family': 9600, 'requested_repetitions': 1000,
            'original_quotas': summary['reference_quotas'], 'pools': {'African': 1000, 'Asian': 2000, 'European': 647},
            'stability_rule': 'q<.05 in >=80% and same finite nonzero original-g sign in >=95% of draws'}
    plan['ladder'] = [summary['reference_quotas'], {'African': 104, 'Asian': 380, 'European': 647},
                      {'African': 208, 'Asian': 760, 'European': 647}, plan['pools']]
    rows = []
    for design, sizes, stable in [('A_original', plan['ladder'][0], [76, 193]),
                                   ('C_stage_1', plan['ladder'][1], [110, 220]),
                                   ('C_stage_2', plan['ladder'][2], [100, 210]),
                                   ('C_stage_3', plan['ladder'][3], [None, None])]:
        for race, count, total in zip(['African', 'Asian'], stable, [135, 252]):
            rows.append({'design': design, 'canonical_race': race, 'reference_n': sizes[race],
                         'repetitions': 1000, 'reference_selection_varies': count is not None,
                         'sampling_scheme': 'original_without_replacement_subset' if design == 'A_original'
                                            else 'without_replacement_subset',
                         'original_discoveries': total, 'stable_original': count,
                         'stable_original_percent': None if count is None else 100*count/total})
    return v.pd.DataFrame(rows), plan


def test_larger_subset_selected_by_size_not_outcome(summary, comparison_source):
    result = v.comparison_values(summary, *comparison_source)
    assert result['design'] == 'C_stage_2'
    assert result['reference_n'] == {'African': 208, 'Asian': 760, 'European': 647}
    assert (result['stable'], result['total']) == (310, 387)
    assert result['percent'] == pytest.approx(100*310/387)


@pytest.mark.parametrize('row,field,value', [
    (0, 'stable_original', 77),
    (0, 'repetitions', 999),
    (0, 'stable_original_percent', 99),
    (4, 'reference_n', 209),
    (4, 'original_discoveries', 134),
    (4, 'stable_original', 100.5),
    (4, 'stable_original_percent', 99),
    (4, 'reference_selection_varies', False),
    (4, 'sampling_scheme', 'with_replacement_bootstrap'),
])
def test_incompatible_reference_comparison_fails(summary, comparison_source, row, field, value):
    table, plan = comparison_source
    table.loc[row, field] = value
    with pytest.raises(ValueError):
        v.comparison_values(summary, table, plan)


def test_duplicate_reference_rows_fail(summary, comparison_source):
    table, plan = comparison_source
    with pytest.raises(ValueError, match='Duplicate'):
        v.comparison_values(summary, v.pd.concat([table, table.iloc[[0]]]), plan)


def test_no_larger_variable_subset_fails(summary, comparison_source):
    table, plan = comparison_source
    plan['ladder'] = [plan['original_quotas'], plan['pools']]
    with pytest.raises(ValueError, match='No larger resampleable'):
        v.comparison_values(summary, table, plan)


def test_wrong_diagnostic_baseline_fails(summary, comparison_source):
    table, plan = comparison_source
    plan['baseline_id'] = 'OLD_POOLED_STUDY'
    with pytest.raises(ValueError, match='baseline/family mismatch'):
        v.comparison_values(summary, table, plan)


def test_released_aggregates_reproduce_plotted_values(summary):
    result = v.load_comparison(summary, source='released')
    plot = v.pd.read_csv(v.RELEASE_INPUT.parent/'panel_C_reference_comparison_values.csv')
    assert len(plot) == 4
    for race in ['African', 'Asian']:
        row = plot.loc[plot.canonical_race.eq(race) & plot.design.eq('Larger subset')].iloc[0]
        assert result['groups'][race] == {'stable': row.stable_original, 'total': row.original_discoveries}
        assert result['reference_n'][race] == row.reference_n
    assert result['input_sha256'][str(v.RELEASE_INPUT.relative_to(v.ROOT))] == v.EXPECTED_RELEASE_INPUT_SHA256


def test_changed_released_aggregate_fails(summary, tmp_path, monkeypatch):
    path = tmp_path/'reference_comparison_inputs.json'
    path.write_text(v.RELEASE_INPUT.read_text()+' ')
    monkeypatch.setattr(v, 'RELEASE_INPUT', path)
    with pytest.raises(ValueError, match='Released reference plotting input changed'):
        v.load_comparison(summary, source='released')


def test_panel_c_preserves_original_and_adds_larger_varied_groups(summary, comparison_source):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.5,4.5))
    values = v.values(summary)
    for race in ['African', 'Asian']:
        values[race.lower() + '_reference_n'] = summary['reference_quotas'][race]
    v.panel_c(ax,values,v.comparison_values(summary, *comparison_source))
    assert len(ax.patches)==4
    assert [patch.get_width() for patch in ax.patches]==pytest.approx([100*76/135,100*100/135,100*193/252,100*210/252])
    assert any('Original: 76/135' in t.get_text() for t in ax.texts)
    assert any('Larger subset: 100/135' in t.get_text() for t in ax.texts)
    assert not any('European recorded category' in t.get_text() for t in ax.texts)
    plt.close(fig)


def test_panel_b_overlap_is_not_mc_count(summary):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.5,3.4))
    v.panel_b(ax,v.values(summary))
    assert [patch.get_width() for patch in ax.patches]==[80,23,58,34]
    assert any('0 pass simultaneous MC upper bounds' in t.get_text() for t in ax.texts)
    plt.close(fig)
