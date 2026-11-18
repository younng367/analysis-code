"""Read-only source extraction and reproducible revision analyses.

Run using tbi_omicverse Python. No source file is overwritten.
Current-cohort animal data are pending and excluded from this analysis.
"""
from pathlib import Path
import io, json, zipfile, hashlib, math, itertools
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats, sparse
from scipy.io import mmwrite
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
BIO=ROOT.parents[1]/'TBI_es_补'
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
sources=[]
def source(p):
    p=Path(p); sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return p
# Pseudobulk is based on raw counts summed within each library, not nuclei.
f=source(BIO/'Micro_ECM.h5ad');a=ad.read_h5ad(f)
counts=a.layers['counts'];assert np.all(counts.data==np.floor(counts.data))
samples=sorted(a.obs['orig.ident'].unique());agg=[];meta=[]
for s in samples:
    m=np.asarray(a.obs['orig.ident']==s);agg.append(np.asarray(counts[m].sum(axis=0)).ravel())
    meta.append({'sample':s,'group':str(a.obs.loc[m,'group'].iloc[0]),'n_nuclei':int(m.sum())})
pd.DataFrame(np.array(agg).T,index=a.var_names,columns=samples).to_csv(OUT/'microglia_pseudobulk_counts.csv.gz',compression='gzip')
pd.DataFrame(meta).to_csv(OUT/'microglia_library_metadata.csv',index=False)
obs=a.obs[['orig.ident','group','cell_type','ECMdis_group','MM_ECMdis_score']].copy()
obs.columns=['sample','timepoint','state','original_score_group','ecm_score']
assert np.array_equal(np.where(obs.ecm_score>=obs.ecm_score.median(),'High','Low'),obs.original_score_group.astype(str))
mm=pd.read_excel(source(BIO/'BULK/bulk/MECHANIC GENE.xlsx'),sheet_name='Mechanical Microenvironment')
gs=[str(g).strip().capitalize() for g in mm['Extracellular matrix disassembly'].dropna()]
gs=[g for g in gs if g in a.var_names]
pd.DataFrame({'gene':gs,'excluded_for_sensitivity':[g in ['Lcp1','Coro1a'] for g in gs]}).to_csv(OUT/'ECM_disassembly_gene_panel.csv',index=False)
sc.tl.score_genes(a,gene_list=gs,score_name='revision_score_reproduction',ctrl_size=50,random_state=0)
sc.tl.score_genes(a,gene_list=[g for g in gs if g not in ['Lcp1','Coro1a']],gene_pool=[g for g in a.var_names if g not in ['Lcp1','Coro1a']],score_name='revision_score_excluding_targets',ctrl_size=50,random_state=0)
obs['ecm_excluding_targets']=a.obs['revision_score_excluding_targets'].to_numpy()
pd.DataFrame([{'scanpy_version':sc.__version__,'n_matched_genes':len(gs),'max_absolute_difference_from_stored_score':float(np.max(np.abs(a.obs.revision_score_reproduction-obs.ecm_score))),'targets_in_score':';'.join(g for g in gs if g in ['Lcp1','Coro1a'])}]).to_csv(OUT/'ECM_score_reproduction_check.csv',index=False)
for gene in ['Lcp1','Coro1a']:
    v=a[:,gene].X;obs[gene]=np.asarray(v.toarray() if sparse.issparse(v) else v).ravel()
obs.to_csv(OUT/'microglia_scores_and_expression.csv.gz',compression='gzip')
obs.groupby(['sample','timepoint','state'],observed=True).agg(n=('ecm_score','size'),ecm_mean=('ecm_score','mean'),ecm_median=('ecm_score','median'),Lcp1_mean=('Lcp1','mean'),Coro1a_mean=('Coro1a','mean')).reset_index().to_csv(OUT/'microglia_sample_state_summary.csv',index=False)
models=[]
for gene in ['Lcp1','Coro1a']:
    for score in ['ecm_score','ecm_excluding_targets']:
        model=smf.ols(f'{gene} ~ {score} + C(sample) + C(state)',data=obs).fit(cov_type='cluster',cov_kwds={'groups':obs['sample']},use_t=True)
        ci=model.conf_int().loc[score];models.append({'gene':gene,'score':score,'model':'log-normalized expression ~ continuous ECM score + library fixed effects + state fixed effects; library-clustered SE with t reference','n_libraries':len(samples),'coefficient':model.params[score],'SE':model.bse[score],'ci95_low':ci.iloc[0],'ci95_high':ci.iloc[1],'p':model.pvalues[score],'interpretation':'exploratory association; 12 clusters; not CellChat robustness'})
pd.DataFrame(models).to_csv(OUT/'ECM_continuous_models.csv',index=False)
thresholds=[]
for rule in ['global_median','within_library_median','global_tertiles','global_quartiles']:
    v=obs.ecm_score
    if rule=='within_library_median':lo=hi=obs.groupby('sample',observed=True).ecm_score.transform('median')
    elif rule=='global_median':lo=hi=v.median()
    else:
        q=1/3 if rule=='global_tertiles' else .25;lo=v.quantile(q);hi=v.quantile(1-q)
    cat=np.where(v>=hi,'High',np.where(v<lo,'Low','Middle'))
    t=obs.assign(score_group=cat).groupby(['sample','timepoint','state','score_group'],observed=True).size().reset_index(name='n')
    t['rule']=rule;thresholds.append(t)
pd.concat(thresholds).to_csv(OUT/'ECM_threshold_composition.csv',index=False)

# Reproduce all network scores, retaining rank ties and isolated candidates.
b=BIO/'BULK';edge=pd.read_csv(source(b/'String/string612.tsv'),sep='\t');cand=pd.read_csv(source(b/'common_genes.csv')).iloc[:,0].tolist()
name_map={s.upper():s for s in set(edge.iloc[:,0])|set(edge.iloc[:,1])}
G=nx.Graph();G.add_nodes_from([name_map.get(c,c) for c in cand]);G.add_edges_from(zip(edge.iloc[:,0],edge.iloc[:,1]))
rank=[]
for gene in G:
    sub=G.subgraph(list(G.neighbors(gene)));comps=list(nx.connected_components(sub));largest=max(comps,key=len) if comps else []
    n=len(largest);e=sub.subgraph(largest).number_of_edges()
    cliques=[c for c in nx.find_cliques(G) if gene in c and len(c)>1]
    rank.append({'gene':gene,'degree':G.degree(gene),'MNC':n,'DMNC':e/(n**1.7) if n else 0,'MCC':sum(math.factorial(len(c)-1) for c in cliques)})
r=pd.DataFrame(rank)
for metric in ['DMNC','MCC','MNC']:r[metric+'_rank']=r[metric].rank(ascending=False,method='min').astype(int)
r['top3_all_including_ties']=(r[['DMNC_rank','MCC_rank','MNC_rank']]<=3).all(axis=1)
r.to_csv(OUT/'original_network_full_rankings.csv',index=False)
rob=[]
for k in range(2,9):
    keep=r[(r[['DMNC_rank','MCC_rank','MNC_rank']]<=k).all(axis=1)]
    rob.append({'rank_cutoff':k,'rule':'minimum rank <= cutoff in all three metrics, retaining ties','genes':';'.join(keep.gene),'n':len(keep)})
pd.DataFrame(rob).to_csv(OUT/'original_network_rank_sensitivity.csv',index=False)
(OUT/'source_hashes_01.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
print('Pseudobulk libraries',len(samples),'top3 all',r.loc[r.top3_all_including_ties,'gene'].tolist())
