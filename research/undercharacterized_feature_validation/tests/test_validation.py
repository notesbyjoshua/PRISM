import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

PATH = Path(__file__).resolve().parents[1] / 'validate_candidates.py'
spec = importlib.util.spec_from_file_location('candidate_validation', PATH)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_median_ties_and_percentile_orientation():
    np.testing.assert_allclose(v.percentiles(np.array([0, 1, 2, 3, 4]), np.array([1, 2, 2, 3])), [0, 12.5, 50, 87.5, 100])


def test_effect_and_u_orientation():
    x, y = np.arange(10)+10., np.arange(20)*.1
    assert v.core.effect(x,y)>0 and v.core.effect(y,x)<0
    assert mannwhitneyu(x,y).statistic/(len(x)*len(y)) == 1
    assert mannwhitneyu(x,y).pvalue == mannwhitneyu(y,x).pvalue


def test_every_case_omitted_once_fixed_reference():
    x,y=np.array([1.,2.,3.,4.,12.]),np.array([0.,1.,2.,3.,4.,5.])
    result=v.loo(x,y)
    assert len(result)==len(x)
    for i,row in enumerate(result.itertuples()):
        expected=np.delete(x,i)
        assert row.g==v.core.effect(expected,y)
        assert row.p_raw==mannwhitneyu(expected,y,method='auto',alternative='two-sided').pvalue
    assert not any('q' in col for col in result)


def test_duplicates_rejected():
    with pytest.raises(ValueError,match='Distinct recorded patient'):
        v.check_patients(pd.DataFrame({'patient_id':['a','a']}))


def test_missing_components_are_not_inferred():
    assert v.decomposition(None,None,np.array([1.]))['status']=='not_evaluable'


def test_ratio_same_numerator_different_denominator():
    numerator=np.array([10.,10.])
    denominator=np.array([2.,4.])
    ratio=np.array([5.,2.5])
    assert v.decomposition(numerator,denominator,ratio)['status']=='formula_verified'
    assert not v.decomposition(numerator,denominator,ratio)['absolute_anatomy_established']
    with pytest.raises(AssertionError):
        v.decomposition(numerator,denominator,np.array([5.,5.]))
    with pytest.raises(ValueError):
        v.decomposition(numerator,np.array([0.,4.]),ratio)


def test_missing_candidate_input_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        v.locked_candidates(tmp_path/'missing.md')
    bad=tmp_path/'wrong.md'
    bad.write_text('No candidates')
    with pytest.raises(ValueError,match='Exactly ten'):
        v.locked_candidates(bad)


def test_private_output_guard(tmp_path):
    with pytest.raises(ValueError,match='inside local_review'):
        v.core.baseline.private_output(tmp_path)


def test_full_family_bh_not_ten_test_adjustment():
    ten=np.repeat(.0001,10)
    full=np.r_[ten,np.ones(9590)]
    assert np.all(v.core.bh(ten)<v.core.bh(full)[:10])


def test_plan_fixed_and_european_na():
    import json
    plan=json.loads((PATH.parent/'analysis_plan.json').read_text())
    assert plan['version']=='1.0.0' and plan['case_bootstrap_B']==2000
    assert 'no variation' in plan['comparator_stability_score']['NA']
    assert 'no LOO BH' in plan['loo']
