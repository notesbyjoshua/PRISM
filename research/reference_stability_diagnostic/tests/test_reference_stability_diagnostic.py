import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import mannwhitneyu

p=Path(__file__).resolve().parents[1]/'engine.py'
spec=importlib.util.spec_from_file_location('diagnostic_engine',p)
e=importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)


@pytest.mark.parametrize('n,m',[(3,9),(5,20),(8,31),(6,647),(5,1174),(8,3228)])
def test_exact_probability_matches_public_scipy(n,m):
    for shift in [0,m//3,m//2,m-1]:
        x=np.arange(n)+shift
        y=np.setdiff1d(np.arange(n+m),x)
        result=mannwhitneyu(x,y,method='exact',alternative='two-sided')
        u=result.statistic
        p=e.exact_lookup(n,m)[int(min(u,n*m-u))]
        np.testing.assert_allclose(p,result.pvalue,rtol=5e-10,atol=1e-14)


@pytest.mark.parametrize('x,y',[(np.array([1,2,2,3,4.]),np.array([1,1,2,2,3,4,4,5.])),
                               (np.arange(8.),np.arange(8.,20.)),(np.ones(6),np.ones(21))])
def test_mw_ties_no_effect_and_shifted_cases(x,y):
    result=mannwhitneyu(x,y,alternative='two-sided',method='auto')
    _,counts=np.unique(np.r_[x,y],return_counts=True)
    tie=np.sum(counts**3-counts)
    p=e.mw_p(np.array([[result.statistic]]),np.array([[len(x)]]),len(y),np.array([[tie]]))[0,0]
    np.testing.assert_allclose(p,result.pvalue,rtol=1e-12,atol=1e-15)


def test_near_threshold_sign_reversal_and_nan_stability():
    g=np.array([[1,1,1,np.nan,0],[1,1,-1,1,0],[1,1,1,1,0],[1,1,1,1,0],[1,1,1,1,0.]])
    q=np.array([[.01,.05,.01,.01,.01]]*5)
    frequency,sign,stable=e.stability(g,q,np.array([1,1,1,1,0.]))
    assert stable.tolist()==[True,False,False,False,False]
    assert frequency[1]==0 and sign[2]==.8 and sign[3]==.8 and sign[4]==0


def test_ladder_full_pool_is_one_subset():
    quotas={'African':5,'Asian':9,'European':20}
    pools={'African':23,'Asian':40,'European':20}
    stages=e.ladder(quotas,pools)
    assert stages==[{'African':5,'Asian':9,'European':20},{'African':10,'Asian':18,'European':20},
                    {'African':20,'Asian':36,'European':20},pools]
    assert len({e.digest(np.ones(20)) for _ in range(10)})==1


def test_weighted_bank_matches_scipy_for_duplicates_and_shared_values():
    import pandas as pd
    from research.race_stratified_sensitivity import core
    rng=np.random.default_rng(72)
    features=['fabricated_a','fabricated_b']
    patient_rows=[];reference_rows=[];cell_rows=[]
    for race in core.RACES:
        y=np.c_[np.round(rng.normal(size=30)),rng.normal(size=30)]
        reference_rows.extend({'file':f'fabricated_{race}_{i}','canonical_race':race,**dict(zip(features,row))} for i,row in enumerate(y))
        for syndrome,n in [('fake_A',5),('fake_B',9)]:
            cell_rows.append({'disease':syndrome,'canonical_race':race,'eligible':True,'n_current_gmdb':n})
            x=y[:n].copy()
            patient_rows.extend({'disease':syndrome,'canonical_race':race,'patient_id':f'{race}_{syndrome}_{i}',
                                 **dict(zip(features,row))} for i,row in enumerate(x))
    counts=pd.DataFrame(cell_rows)
    cells=counts.sort_values(['disease','canonical_race'])
    primary=pd.DataFrame([{'disease':row.disease,'canonical_race':row.canonical_race,'feature':f}
                          for row in cells.itertuples() for f in features])
    loaded=(pd.DataFrame(patient_rows),pd.DataFrame(reference_rows),features,None,{},counts,primary,{},[])
    bank=e.Bank(loaded,core)
    for indices in [np.arange(20),np.tile(np.arange(15),2),np.arange(30)]:
        draw={r:indices for r in core.RACES}
        p,g,_,_=bank.evaluate(draw)
        for race,data in bank.races.items():
            for j,x in enumerate(data['cases']):
                _,expected_p,expected_g=core.tests_direct(x,data['y'][indices])
                where=data['where'].reshape(len(data['cases']),-1)[j]
                np.testing.assert_allclose(p[where],expected_p,rtol=1e-12,atol=1e-15)
                np.testing.assert_allclose(g[where],expected_g,rtol=1e-12,atol=1e-15)


def test_atomic_checkpoint_resume_signature_and_corruption(tmp_path):
    from research.race_stratified_sensitivity import core
    cp=core.Checkpoints(tmp_path,'fake','known_signature',False)
    cp.append(2,g=np.ones((2,12)),q=np.ones((2,12)))
    restored=core.Checkpoints(tmp_path,'fake','known_signature',True)
    assert restored.state['completed']==2
    with pytest.raises(ValueError,match='signature changed'):
        core.Checkpoints(tmp_path,'fake','different_signature',True)
    file=cp.directory/cp.state['chunks'][0]['file']
    file.write_bytes(b'corrupted fabricated fixture')
    with pytest.raises(ValueError,match='corruption'):
        core.Checkpoints(tmp_path,'fake','known_signature',True)


def test_no_effect_and_reversed_effect_have_expected_signs():
    from research.race_stratified_sensitivity import core
    x=np.arange(10,dtype=float)
    assert core.effect(x,x)==0
    assert core.effect(x+10,x)>0 and core.effect(x-10,x)<0
    frequency,sign,stable=e.stability(np.array([[1.],[1.]]),np.array([[np.nan],[np.nan]]),np.array([1.]))
    assert frequency[0]==0 and not stable[0]
