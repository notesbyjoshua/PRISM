import importlib.util
from pathlib import Path

import numpy as np

path=Path(__file__).resolve().parents[1]/'recommend.py'
spec=importlib.util.spec_from_file_location('recommendations',path)
v=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_category_baselines_cancel_for_equal_case_effects():
    value,se,df,stat,p,status=v.contrast([12,10,102,100],[1,2,1,2],[10,20,10,20])
    assert value==0 and p==1 and stat==0 and status=='tested'


def test_four_group_variance_and_satterthwaite_df():
    ns=np.array([6,52,8,190])
    sds=np.array([1.,2.,3.,4.])
    value,se,df,stat,p,status=v.contrast([5,2,2,1],sds,ns)
    variances=sds**2/ns
    np.testing.assert_allclose(se**2,variances.sum())
    np.testing.assert_allclose(df,variances.sum()**2/np.sum(variances**2/(ns-1)))
    assert value==2 and stat>0 and 0<p<1


def test_swapping_categories_reverses_sign_only():
    a=v.contrast([6,1,3,2],[1,2,3,4],[10,20,30,40])
    b=v.contrast([3,2,6,1],[3,4,1,2],[30,40,10,20])
    np.testing.assert_allclose(a[0],-b[0])
    np.testing.assert_allclose(a[4],b[4])


def test_zero_variance_has_no_significance_claim():
    result=v.contrast([1,0,3,0],[0,0,0,0],[5,20,5,20])
    assert result[-1]=='undefined_zero_variance' and result[4]==1
