import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from research.race_stratified_pitch_sensitivity import design as d
from research.race_stratified_pitch_sensitivity.reporting import figures
from research.race_stratified_pitch_sensitivity.validate import independent_tests


@pytest.fixture
def fabricated():
    rng = np.random.default_rng(332)
    features = [f'fabricated_feature_{i}' for i in range(120)]
    cells = pd.DataFrame({'disease': [f'Fabricated_{i}' for i in range(80)],
                          'canonical_race': ['African','Asian']+['European']*78,
                          'n_current_gmdb': [10,10]+[5]*78})
    pools = {'African': (200,100,100), 'Asian': (200,100,100), 'European': (647,323,324)}
    quotas = {'African':20,'Asian':20,'European':647}
    support = d.support_table(cells,pools,quotas)
    groups = {}
    for race, (n,_,_) in pools.items():
        pitch = np.linspace(-14,14,n)
        keys = np.asarray([f'fabricated_key_{i:04d}' for i in range(n)])
        low, high = d.pitch_halves(pitch,keys,'degrees')
        groups[race] = {'values': rng.normal(size=(n,120)), 'pose': np.c_[pitch,rng.normal(size=(n,2))],
                        'sex': np.where(np.arange(n)%2,'Female','Male'), 'low':low,'high':high}
    return {'features':features,'cells':cells,'groups':groups,'support':support,
            'supported':np.repeat(support.supported.to_numpy(bool),120),'quotas':quotas}


def test_pitch_halves_use_magnitude_lexical_ties_and_odd_high_assignment():
    pitch, keys = np.array([-3,1,-1,0,2]), np.array(['c','b','a','e','d'])
    low, high = d.pitch_halves(pitch,keys,'degrees')
    np.testing.assert_array_equal(low,[3,2])
    np.testing.assert_array_equal(high,[1,4,0])
    order = np.array([3,1,4,0,2])
    shuffled_low, shuffled_high = d.pitch_halves(pitch[order],keys[order],'degrees')
    assert set(keys[order][shuffled_low]) == set(keys[low])
    assert set(keys[order][shuffled_high]) == set(keys[high])
    assert not set(low)&set(high) and set(low)|set(high) == set(range(5))


@pytest.mark.parametrize('pitch,units', [([0,np.nan],'degrees'),([0,np.inf],'degrees'),([0,1],'unknown')])
def test_nonfinite_or_unverified_pitch_is_blocked(pitch,units):
    with pytest.raises(d.Blocked,match='BLOCKED'):
        d.pitch_halves(pitch,['a','b'],units)


@pytest.mark.parametrize('fault',['missing_pitch','nonfinite_pitch','unverified_units','old_QC_violation'])
def test_pose_contract_has_documented_blocked_audit(fault):
    frame = pd.DataFrame({'canonical_race':['African','Asian','European'],
                          'pose_pitch':[1.,2.,3.],'pose_yaw':[0.,0.,0.],'pose_roll':[0.,0.,0.]})
    units = 'degrees'
    if fault == 'missing_pitch':
        frame = frame.drop(columns='pose_pitch')
    elif fault == 'nonfinite_pitch':
        frame.loc[0,'pose_pitch'] = np.nan
    elif fault == 'unverified_units':
        units = 'radians_unknown'
    else:
        frame.loc[0,'pose_pitch'] = 16
    with pytest.raises(d.Blocked) as error:
        d.require_pose(frame,{'pitch':15,'yaw':15,'roll':10},units)
    assert str(error.value).startswith('BLOCKED:') and len(error.value.audit) == 9


@pytest.mark.parametrize('arm',d.ARMS)
@pytest.mark.parametrize('capacity,sizes,reference_n,expected',[(200,[5,9],20,'shared_disjoint'),(60,[22,22],20,'per_pair_disjoint_fallback')])
def test_exact_disjoint_sampling_and_reproducibility(arm,capacity,sizes,reference_n,expected):
    low, high = np.arange(capacity//2), np.arange(capacity//2,capacity)
    first = d.allocations(capacity,low,high,sizes,reference_n,d.rng_for(418,12,0,d.ARMS.index(arm)),arm)
    second = d.allocations(capacity,low,high,sizes,reference_n,d.rng_for(418,12,0,d.ARMS.index(arm)),arm)
    assert first[0] == second[0] == expected
    for size,pair,repeat in zip(sizes,first[1],second[1]):
        case, reference = pair
        assert len(case) == len(set(case)) == size and len(reference) == len(set(reference)) == reference_n
        assert not set(case)&set(reference)
        np.testing.assert_array_equal(case,repeat[0])
        np.testing.assert_array_equal(reference,repeat[1])
        if arm == 'pitch_imbalanced':
            assert set(case).issubset(set(high)) and set(reference).issubset(set(low))


def test_substreams_do_not_depend_on_execution_order_or_chunk_size():
    forward = {i:d.rng_for(123,i,1,0).integers(100000,size=10).tolist() for i in range(30)}
    reverse = {i:d.rng_for(123,i,1,0).integers(100000,size=10).tolist() for i in reversed(range(30))}
    chunked = {i:d.rng_for(123,i,1,0).integers(100000,size=10).tolist() for start in [0,7,14,21,28] for i in range(start,min(start+7,30))}
    assert forward == reverse == chunked
    assert forward[0] != d.rng_for(123,0,1,1).integers(100000,size=10).tolist()


def test_common_support_and_new_split_shortages(fabricated):
    support = fabricated['support']
    assert support.supported.sum() == 2 and np.array_equal(support.supported,support.nominal_random_supported)
    assert support.loc[support.canonical_race.eq('European'),'support_reason'].str.contains('random_insufficient_disjoint_pool').all()
    pools = {'African':(200,100,100),'Asian':(60,30,30),'European':(647,323,324)}
    quotas = {'African':20,'Asian':40,'European':647}
    changed = d.support_table(fabricated['cells'],pools,quotas)
    assert changed.loc[1,'nominal_random_supported'] and changed.loc[1,'newly_unsupported_after_split']
    assert not changed.loc[1,'supported'] and 'low_pitch_insufficient' in changed.loc[1,'support_reason']


def test_full_family_keeps_unsupported_p_one_and_never_discovers(fabricated):
    supported = fabricated['supported']
    p,g,q,found = d.full_family(np.zeros(supported.sum()),np.ones(supported.sum()),supported,fabricated['features'])
    assert len(p) == len(q) == 9600
    assert (p[~supported] == 1).all() and (q[~supported] == 1).all() and not found[~supported].any()
    assert np.isnan(g[~supported]).all()
    np.testing.assert_allclose(q,multipletests(p,method='fdr_bh')[1])
    with pytest.raises(ValueError):
        d.full_family(p[supported],g[supported],supported[:-1],fabricated['features'])


def test_feature_schema_rejects_swapped_columns(fabricated):
    features = fabricated['features']
    frame = pd.DataFrame(np.ones((3,120)),columns=features)
    assert d.locked.baseline.validate_features(frame,features).shape == (3,120)
    with pytest.raises(ValueError,match='Feature order'):
        d.locked.baseline.validate_features(frame[features[::-1]],features)


@pytest.mark.parametrize('n,m',[(5,20),(8,25),(9,21),(12,28)])
@pytest.mark.parametrize('ties',[False,True])
def test_independent_scipy_and_exact_gamma_effect(n,m,ties):
    rng = np.random.default_rng(123+n)
    x,y = rng.normal(size=(n,5)),rng.normal(size=(m,5))
    if ties:
        x[:,0], y[:,0] = np.round(x[:,0]),np.round(y[:,0])
        x[:,1], y[:,1] = 1,1
    actual = d.locked.tests_direct(x,y)
    expected = independent_tests(x,y)
    for first,second in zip(actual,expected):
        np.testing.assert_allclose(first,second,rtol=1e-10,atol=1e-12,equal_nan=True)
    if ties:
        assert actual[1][1] == 1 and np.isnan(actual[2][1])


def test_designed_pitch_feature_contrast_and_independent_null(fabricated):
    for group in fabricated['groups'].values():
        group['values'] = np.repeat(np.abs(group['pose'][:,0,None]),120,axis=1)
    discoveries = {arm:[] for arm in d.ARMS}
    for rep in range(12):
        for arm_index,arm in enumerate(d.ARMS):
            p,g,pose,_ = d.evaluate(fabricated,882,rep,arm_index)
            discoveries[arm].append(d.full_family(p,g,fabricated['supported'],fabricated['features'])[3].sum())
            if arm_index == 1:
                assert (pose[:,12:14] > 0).all() and (pose[:,14] == 0).all()
    assert np.mean(discoveries['pitch_imbalanced']) > np.mean(discoveries['random'])+100
    rng = np.random.default_rng(78)
    for group in fabricated['groups'].values():
        group['values'] = rng.normal(size=group['values'].shape)
    independent = []
    for rep in range(20):
        p,g,_,_ = d.evaluate(fabricated,882,rep,1)
        independent.append(d.full_family(p,g,fabricated['supported'],fabricated['features'])[3].sum())
    assert np.mean(independent) < 5


def test_constant_pitch_cannot_be_reported_as_achieved_imbalance(fabricated):
    for group in fabricated['groups'].values():
        group['pose'][:,0] = 0
    with pytest.raises(d.Blocked,match='Achieved pitch separation'):
        d.evaluate(fabricated,7,0,1)


def test_checkpoint_requires_matching_signature_and_checks_hashes(tmp_path):
    checkpoint = d.locked.Checkpoints(tmp_path,'fabricated','one',False)
    checkpoint.append(2,p=np.ones((2,2,240)))
    with pytest.raises(ValueError,match='requires --resume'):
        d.locked.Checkpoints(tmp_path,'fabricated','one',False)
    with pytest.raises(ValueError,match='signature changed'):
        d.locked.Checkpoints(tmp_path,'fabricated','two',True)
    assert d.locked.Checkpoints(tmp_path,'fabricated','one',True).state['completed'] == 2
    path=checkpoint.directory/checkpoint.state['chunks'][0]['file']
    path.write_bytes(path.read_bytes()+b'corruption')
    with pytest.raises(ValueError,match='corruption'):
        d.locked.Checkpoints(tmp_path,'fabricated','one',True)


def test_figures_render_and_mark_european_unavailable(tmp_path):
    rng = np.random.default_rng(49)
    reps = pd.DataFrame({'arm':np.repeat(d.ARMS,10),'discoveries':np.r_[rng.poisson(.1,10),rng.poisson(100,10)]})
    categories, pose = [],[]
    for arm in d.ARMS:
        for race in d.RACES:
            supported = race != 'European'
            values = reps.loc[reps.arm.eq(arm),'discoveries'].to_numpy()
            categories.append({'arm':arm,'canonical_race':race,'supported_tests':120 if supported else 0,
                               'runs':10 if supported else np.nan,'runs_with_any':np.sum(values>0) if supported else np.nan,
                               'fraction_with_any':np.mean(values>0) if supported else np.nan})
            pose.append({'arm':arm,'canonical_race':race,'metric':'mean_abs_pitch_gap',
                         'mean':(0.02 if arm == 'random' else 7.0) if supported else np.nan})
    figures(tmp_path,reps,pd.DataFrame(categories),pd.DataFrame(pose),2,240,647)
    for stem in ['pitch_placebo_distribution','pitch_imbalance_diagnostics']:
        assert all((tmp_path/(stem+'.'+ext)).is_file() for ext in ['png','pdf','svg'])
    svg=(tmp_path/'pitch_imbalance_diagnostics.svg').read_text()
    assert 'N/A' in svg and '647' in svg and '<image ' not in svg
