from pathlib import Path
import json,itertools,gzip,concurrent.futures
import numpy as np,pandas as pd,requests
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
s=pd.read_csv(OUT/'bulk_immune_ssGSEA_scores.csv',index_col=0)
e=pd.read_csv(OUT/'bulk_normalized_counts.csv',index_col=0)
assert list(s.columns)==list(e.columns)
rows=[];points=[]
perms=np.array(list(itertools.permutations(range(6))))
for gene in ['Lcp1','Coro1a']:
    x=e.loc[gene].values;xr=rankdata(x);xr=(xr-xr.mean())/np.linalg.norm(xr-xr.mean())
    for signature,y in s.iterrows():
        yr=rankdata(y.values);yr=(yr-yr.mean())/np.linalg.norm(yr-yr.mean());rho=float(xr@yr)
        p=float(np.mean(np.abs(yr[perms]@xr)>=abs(rho)-1e-12))
        rows.append({'gene':gene,'signature':signature,'n_samples':6,'rho':rho,'p_exact_permutation':p})
        for sample,v,score in zip(s.columns,x,y.values):points.append({'gene':gene,'signature':signature,'sample':sample,'group':'TBI' if sample.startswith('TBI') else 'Control','normalized_count':v,'ssGSEA_score':score})
r=pd.DataFrame(rows);r['BH_56_tests']=multipletests(r.p_exact_permutation,method='fdr_bh')[1]
for gene in r.gene.unique():
    idx=r.gene==gene;r.loc[idx,'BH_within_gene_28_tests']=multipletests(r.loc[idx,'p_exact_permutation'],method='fdr_bh')[1]
r.to_csv(OUT/'bulk_correlations_exact.csv',index=False);pd.DataFrame(points).to_csv(OUT/'bulk_correlation_points.csv',index=False)

def get_geo(acc):
    prefix=acc[:-3]+'nnn';url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/soft/{acc}_family.soft.gz'
    try:
        res=requests.get(url,timeout=90);res.raise_for_status();data=gzip.decompress(res.content).decode()
        (ROOT/'work/intake'/f'{acc}_family.soft').write_text(data)
        kept=[l for l in data.splitlines() if l.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics_ch1','!Sample_description','!Sample_source_name_ch1','!Sample_growth_protocol_ch1','!Sample_treatment_protocol_ch1','!Series_title','!Series_pubmed_id'))]
        (ROOT/'work/intake'/f'{acc}_metadata_excerpt.txt').write_text('\n'.join(kept))
        return {'accession':acc,'url':url,'status':'retrieved','lines':len(kept)}
    except Exception as ex:return {'accession':acc,'url':url,'error':str(ex)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:meta=list(pool.map(get_geo,['GSE79441','GSE277487','GSE236171']))
(OUT/'GEO_source_metadata.json').write_text(json.dumps(meta,indent=2))
print('Correlations',len(r),'BH significant',int((r.BH_56_tests<.05).sum()));print(meta)
