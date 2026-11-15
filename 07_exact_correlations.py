"""Generic exact Spearman permutation comparisons; no study data or file paths."""
from itertools import permutations
from math import factorial

def exact_spearman(x, y, max_permutations=40320):
    import numpy as np
    from scipy.stats import rankdata
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 3:
        raise ValueError('Supply matched one-dimensional sample vectors (n >= 3).')
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError('Missing/non-finite values require an explicit analysis decision.')
    if factorial(len(x)) > max_permutations:
        raise ValueError('Exact enumeration exceeds limit; specify a separate permutation plan.')
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx-rx.mean(), ry-ry.mean()
    denom = np.linalg.norm(rx)*np.linalg.norm(ry)
    if denom == 0:
        raise ValueError('A constant vector has undefined Spearman correlation.')
    rho = float(rx@ry/denom)
    extreme = sum(abs(float(rx@ry[list(p)]/denom)) >= abs(rho)-1e-12
                  for p in permutations(range(len(x))))
    return {'rho':rho, 'exact_two_sided_p':extreme/factorial(len(x)),
            'n_samples':len(x), 'permutations':factorial(len(x))}

def correlate_panel(expression, scores, targets):
    """DataFrames: rows genes/signatures, columns independent samples; BH across all tests."""
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    if not expression.columns.is_unique or not scores.columns.is_unique:
        raise ValueError('Sample names must be unique.')
    if set(expression.columns) != set(scores.columns):
        raise ValueError('Expression and score samples must match.')
    scores = scores.loc[:,expression.columns]
    rows = []
    for target in targets:
        for signature in scores.index:
            rows.append(dict(target=target, signature=signature,
                **exact_spearman(expression.loc[target],scores.loc[signature])))
    result = pd.DataFrame(rows)
    result['BH_p'] = multipletests(result.exact_two_sided_p,method='fdr_bh')[1]
    return result

