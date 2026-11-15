"""Templates for prepared log-normalized AnnData and an exploratory score model."""
def ecm_scores(adata, genes, excluded_targets, control_size, seed):
    import scanpy as sc
    a = adata.copy()
    if not a.var_names.is_unique:
        raise ValueError('Gene identifiers must be unique.')
    panel = [g for g in dict.fromkeys(genes) if g in a.var_names]
    reduced = [g for g in panel if g not in excluded_targets]
    pool = [g for g in a.var_names if g not in excluded_targets]
    if not panel or not reduced:
        raise ValueError('No matched score genes remain.')
    # Explicitly use a.X (prepared log-normalized expression), not raw counts.
    sc.tl.score_genes(a, panel, score_name='score', use_raw=False,
                     ctrl_size=control_size, random_state=seed)
    sc.tl.score_genes(a, reduced, gene_pool=pool, score_name='score_excluding_targets',
                     use_raw=False, ctrl_size=control_size, random_state=seed)
    return a

def continuous_score_model(table):
    """Columns expression, score, sample, state; rows nuclei with library membership."""
    import statsmodels.formula.api as smf
    required = ['expression','score','sample','state']
    if table[required].isna().any().any() or table['sample'].nunique() < 3:
        raise ValueError('Complete values and multiple independent libraries required.')
    fit = smf.ols('expression ~ score + C(sample) + C(state)', data=table).fit(
        cov_type='cluster', cov_kwds={'groups':table['sample']}, use_t=True)
    ci = fit.conf_int().loc['score']
    return {'estimate':fit.params['score'], 'SE':fit.bse['score'],
            'CI_low':ci.iloc[0], 'CI_high':ci.iloc[1], 'p':fit.pvalues['score'],
            'n_libraries':table['sample'].nunique()}

def score_threshold_groups(table):
    """Return labels only; these splits do not validate a CellChat network."""
    import numpy as np
    x = table['score']
    if x.isna().any() or table['sample'].isna().any():
        raise ValueError('Missing scores or sample IDs.')
    within = table.groupby('sample')['score'].transform('median')
    rules = {'global_median':(x.median(),x.median()),
             'within_sample_median':(within,within),
             'global_tertiles':(x.quantile(1/3),x.quantile(2/3)),
             'global_quartiles':(x.quantile(.25),x.quantile(.75))}
    return {name:np.where(x>=hi,'High',np.where(x<lo,'Low','Middle'))
            for name,(lo,hi) in rules.items()}
# A small number of library clusters limits precision; this is an association model.

