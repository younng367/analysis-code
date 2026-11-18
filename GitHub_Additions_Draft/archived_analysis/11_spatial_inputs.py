from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,json,hashlib
from scipy import sparse
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results'
f=ROOT.parents[1]/'TBI_es_补/GSE236171/old/612_crop_st.h5ad';a=ad.read_h5ad(f)
assert a.n_obs==1068 and a.obs.region.value_counts().to_dict()=={'Damaged sites(Core)':590,'Non Damaged sites/Lesion':478}
x=pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']);x['region']=a.obs.region.astype(str)
x.to_csv(out/'spatial_trajectory_input.csv',index_label='spot')
y=pd.DataFrame(a.obsm['spatial'],index=a.obs_names,columns=['x','y']);y['region']=a.obs.region.astype(str)
z=x[['region']].copy()
for gene in ['Lcp1','Coro1a']:
 v=a[:,gene].X;z[gene]=np.asarray(v.toarray() if sparse.issparse(v) else v).ravel();y[gene]=z[gene]
z.to_csv(out/'spatial_target_expression.csv',index_label='spot');y.to_csv(out/'spatial_coordinates.csv',index_label='spot')
(out/'spatial_object_provenance.json').write_text(json.dumps({'file':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'spots':1068,'section_n':1,'timepoint':'3 days post-injury','expression_scale':'saved object X; no additional normalization','original_Slingshot_result':'not recovered; revision refit on saved UMAP'},indent=2))
e=pd.read_csv(out/'ECM_continuous_models.csv');e['BH_two_targets_per_score']=e.groupby('score').p.transform(lambda s:multipletests(s,method='fdr_bh')[1]);e.to_csv(out/'ECM_continuous_models.csv',index=False)
print(x.region.value_counts().to_dict())
