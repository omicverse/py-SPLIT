from pathlib import Path
import numpy as np
import pandas as pd
from anndata import AnnData
import splitst


def test_real_dlpfc_purification_and_reassignment_match_r():
    f = np.load(Path(__file__).parents[1]/'examples/data/split_dlpfc_151673.npz')
    a = AnnData(f['counts'], obs=pd.DataFrame(index=f['cells']), var=pd.DataFrame(index=f['genes']))
    a.layers['counts'] = a.X.copy()
    a.obsm['spatial'] = f['coords']
    w = pd.DataFrame(f['weights'], index=a.obs_names, columns=f['types'])
    r = pd.DataFrame(f['reference'], index=f['types'], columns=a.var_names)
    splitst.purify(a,w,r)
    splitst.spatial_score(a,w,a.obs.first_type,k=12)
    splitst.balance(a,threshold=0.15)
    splitst.reassign_residuals(a,k=12)
    np.testing.assert_allclose(a.layers['split_purified'],f['r_purified'],atol=1e-10,rtol=1e-9)
    np.testing.assert_allclose(a.layers['split_reassigned'].toarray(),f['r_reassigned'],atol=1e-10,rtol=1e-9)
