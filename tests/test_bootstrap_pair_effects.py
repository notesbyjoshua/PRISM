import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'stats_testing'))
from bootstrap_utils import bootstrap_pair_effects


@pytest.mark.parametrize('seed', [1,42,2026])
@pytest.mark.parametrize('n,m', [(10,20),(15,101),(3,4)])
def test_exact_historical_draws_and_effects(seed,n,m):
    rng=np.random.default_rng(100)
    x,y=rng.integers(-3,5,n).astype(float),rng.integers(-2,6,m).astype(float)
    left,right=np.random.default_rng(seed),np.random.default_rng(seed)
    effects,deltas=bootstrap_pair_effects(x,y,207,left)
    expected_g,expected_delta=[],[]
    for _ in range(207):
        a,b=x[right.integers(0,n,n)],y[right.integers(0,m,m)]
        pooled=((n-1)*np.var(a,ddof=1)+(m-1)*np.var(b,ddof=1))/(n+m-2)
        with np.errstate(divide='ignore',invalid='ignore'):
            expected_g.append((1-3/(4*(n+m)-9))*(a.mean()-b.mean())/np.sqrt(pooled))
        expected_delta.append(np.sign(a[:,None]-b[None,:]).sum()/(n*m))
    np.testing.assert_allclose(effects,expected_g,rtol=1e-13,atol=1e-13,equal_nan=True)
    np.testing.assert_array_equal(deltas,expected_delta)
    np.testing.assert_array_equal(left.integers(1000,size=20),right.integers(1000,size=20))


def test_frozen_ci_function_equivalence():
    archived=ROOT/'revision_results/rebuild_120_20261028/local_only/run_20261028/archive/historical_125/stats_testing/statistical_analysis.py'
    old=ast.parse(archived.read_text());new=ast.parse((ROOT/'stats_testing/statistical_analysis.py').read_text())
    names={'cohens_d','hedges_g','cliffs_delta_fast','bootstrap_effect_ci'}
    def functions(tree):
        scope={'np':np,'pd':pd,'bootstrap_pair_effects':bootstrap_pair_effects}
        nodes=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names]
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'<bootstrap-equivalence>','exec'),scope)
        return scope['bootstrap_effect_ci']
    old_fn,new_fn=functions(old),functions(new)
    x=np.array([1.,2.,2.,4.,7.,8.,10.]);y=np.arange(45.)/3
    a,b=np.random.default_rng(2026),np.random.default_rng(2026)
    np.testing.assert_allclose(old_fn(x,y,500,a),new_fn(x,y,500,b),rtol=1e-12,atol=1e-12)
    assert a.integers(100000)==b.integers(100000)
