import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PATH = Path(__file__).resolve().parents[1] / 'generate_facekit_alternative_figures.py'
spec = importlib.util.spec_from_file_location('alternative_figures',PATH)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_fixed_reference_bootstrap_matches_explicit_seeded_loop():
    x = np.array([[1.,2.],[2.,4.],[4.,5.],[5.,7.]])
    y = np.array([[0.,1.],[1.,2.],[2.,2.],[3.,3.],[4.,4.]])
    before = y.copy()
    rng = np.random.default_rng(31)
    indices = rng.integers(0,len(x),size=(200,len(x)))
    expected = np.quantile([v.effect(x[index],y) for index in indices],[.025,.975],axis=0)
    np.testing.assert_allclose(v.bootstrap(x,y,200,31),expected,rtol=1e-12)
    np.testing.assert_array_equal(before,y)
    np.testing.assert_array_equal(v.bootstrap(x,y,200,31),v.bootstrap(x,y,200,31))


def test_bootstrap_failure_never_invents_limits():
    with pytest.raises(ValueError):
        v.bootstrap(np.array([[np.nan],[1.]]),np.array([[1.],[2.]]),10,1)
    with pytest.raises(ValueError):
        v.bootstrap(np.ones((3,1)),np.ones((3,1)),10,1)


def test_selection_never_uses_significance_or_effects():
    rows = [{'disease':d,'canonical_race':r,'case_n':n,'q':.001,'g':200.}
            for d,n in [('a',10),('b',9),('c',8),('d',7),('e',6),('f',5)] for r in ['Asian','European']]
    cells = pd.DataFrame(rows)
    one = v.select(cells).disease.tolist()
    cells['q'] = np.arange(len(cells))/len(cells)
    cells['g'] = -cells['g']
    assert one == v.select(cells).disease.tolist() == ['a','b','c','d','e']


def test_selection_requires_two_categories_and_five_syndromes():
    cells = pd.DataFrame([{'disease':'a','canonical_race':'African','case_n':100}])
    with pytest.raises(ValueError,match='Fewer than five'):
        v.select(cells)


def test_nine_regions_exhaustively_map_to_seven():
    counts = {'forehead':4,'eyebrows':20,'eyes / orbital dimensions':29,'ocular position / iris / gaze':16,
              'nose':17,'philtrum / mouth / lips':15,'chin / jaw':5,'cheeks / malar region':8,
              'global face shape / proportions':6}
    assert set(counts)==set(v.REGION_MAP) and len(set(v.REGION_MAP.values()))==7
    totals = {}
    for region,n in counts.items():
        totals[v.REGION_MAP[region]] = totals.get(v.REGION_MAP[region],0)+n
    assert sum(totals.values())==120 and totals['Eyes / orbital']==45 and totals['Face shape / cheeks']==14


def test_prespecified_features_and_conditional_language():
    assert v.FEATURES == ['eye_fissure_height_mean','nose_tip_width','philtrum_length','chin_width']
    assert 'conditional case-only' in v.FOREST_CAPTION and 'multiplicity' in v.FOREST_CAPTION
    assert 'no region-wise correction' in v.REGION_CAPTION


def test_regions_need_no_raw_records_and_preserve_missing_cells(tmp_path,monkeypatch):
    def forbidden():
        raise AssertionError('Regions must not load patient records')
    monkeypatch.setattr(v.common,'load_sources',forbidden)
    monkeypatch.setattr(v,'save',lambda fig,*args:v.plt.close(fig))
    denominators = dict(zip(v.REGIONS,[4,20,45,17,15,5,14]))
    selected = pd.DataFrame({'disease':[f'Fabricated syndrome {i}' for i in range(5)]})
    rows = []
    for disease in selected.disease:
        for race in ['Asian','European']:
            for region,n in denominators.items():
                rows.extend([{'disease':disease,'canonical_race':race,'display_region':region,
                              'case_n':10,'supported':True,'original_significant':i%2==0} for i in range(n)])
    v.regions(pd.DataFrame(rows),None,selected,denominators,tmp_path)
    out = pd.read_csv(tmp_path/'anatomical_region_chart_data.csv')
    absent = out.loc[out.canonical_race.eq('African')]
    assert absent.status.eq('ineligible').all() and absent.percent_significant.isna().all()
    eligible = out.loc[out.status.eq('eligible')]
    assert len(eligible)==70 and eligible.supported_tests.notna().all()
    assert eligible.groupby(['disease','canonical_race']).planned_features.sum().eq(120).all()
