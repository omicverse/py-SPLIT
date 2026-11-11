from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
import splitst


@pytest.mark.parametrize('threshold', [.05,.15,.3,.5])
@pytest.mark.parametrize('swap', [False,True])
def test_native_r_balance_and_shift(threshold,swap):
    root=Path(__file__).parents[1]
    f=np.load(root/'examples/data/split_dlpfc_151673.npz')
    expected=np.load(root/'tests/data/native_r_balance.npz')
    meta=pd.read_csv(root/'tests/data/native_r_controls.csv',index_col=0)
    for key in ['first_type_neighborhood_agreement','first_type_class_neighborhood_agreement']:
        meta[key]=meta[key].astype(str).str.lower().eq('true')
    a=AnnData(f['counts'].copy(),obs=meta,var=pd.DataFrame(index=f['genes']))
    a.layers['counts']=a.X.copy();a.obsm['spatial']=f['coords']
    w=pd.DataFrame(f['weights'],index=a.obs_names,columns=f['types'])
    ref=pd.DataFrame(f['reference'],index=f['types'],columns=f['genes'])
    splitst.purify(a,w,ref)
    splitst.spatial_score(a,w,a.obs.first_type,k=12)
    splitst.balance(a,threshold=threshold,spot_class_key='spot_class',swap_labels=swap)
    key=f't{threshold}_swap{int(swap)}';cells=expected[key+'_cells']
    actual=pd.DataFrame(a.layers['split_balanced'].toarray(),index=a.obs_names).loc[cells].values
    np.testing.assert_allclose(actual,expected[key+'_counts'],atol=1e-10,rtol=1e-9)
    np.testing.assert_array_equal(a.obs.loc[cells,'split_shift_swap'].values,expected[key+'_swap'])
    np.testing.assert_array_equal(a.obs.loc[cells,'split_balance_status'].values,expected[key+'_status'])
    removed=a.obs.spot_class.eq('reject').values
    assert (a.layers['split_balanced'].toarray()[removed]==0).all()
