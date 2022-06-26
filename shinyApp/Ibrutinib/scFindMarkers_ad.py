import numpy as np
from scipy import stats
from scipy import sparse
from statsmodels.stats.multitest import multipletests
import pandas as pd
import anndata as ad
import gc

def scFindMarkers_py():
    adata = ad.read_h5ad('sc1csr_gexpr.h5ad')
    print(adata.X.shape)
    ident1 = adata[r.Cells1]
    ident2 = adata[r.Cells2]
    ident1 = ident1.X.toarray()
    ident2 = ident2.X.toarray()
    features = adata.var['features']
    del adata
    gc.collect()
    
    N_genes = ident1.shape[1]
    p_vals = np.empty(N_genes)
    p_vals[:] = np.NaN
    for i in range(N_genes):
        _,p_vals[i] = stats.ranksums(ident1[:,i],ident2[:,i])
        
    ident1_pos = ident1 > 0
    pct_1 = np.ravel(ident1_pos.mean(axis=0,dtype=np.float32))
    ident2_pos = ident2 > 0
    pct_2 = np.ravel(ident2_pos.mean(axis=0,dtype=np.float32))

    avg_UMI_1 = np.ravel(np.log2(np.expm1(ident1).mean(axis=0)+1))
    avg_UMI_2 = np.ravel(np.log2(np.expm1(ident2).mean(axis=0)+1))
    
    avg_logFC = avg_UMI_1-avg_UMI_2

    if isinstance(r.inpident_1, list):
        inpident_1 = '_'.join(r.inpident_1)
    else :
        inpident_1 = r.inpident_1

    if isinstance(r.inpident_2, list):
        inpident_2 = '_'.join(r.inpident_2)
    else :
        inpident_2 = r.inpident_2

    markers = pd.DataFrame({'p_val':p_vals,
                            'avg_log2FC':avg_logFC,
                            'pct.1':pct_1,
                            'pct.2':pct_2,
                            #'p_val_adj':p_adj,
                            'avg_log2_UMI.1':avg_UMI_1,
                            'avg_log2_UMI.2':avg_UMI_2,
                            'gene':features,
                            'cluster': ' vs. '.join([inpident_1, inpident_2])})
    markers = markers.dropna()
    _, p_adj, _, _ = multipletests(markers['p_val'],method='fdr_bh')
    markers['p_val_adj'] = p_adj
    markers = markers.sort_values('avg_log2FC', ascending=False)
    markers = markers[['p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj',
                    'avg_log2_UMI.1','avg_log2_UMI.2','gene','cluster']]
    store = pd.HDFStore('tempData/scFindMarkers.h5', 'w', complib=str('zlib'))
    store.put('de', markers, data_columns=markers.columns)
    store.close()

scFindMarkers_py()
