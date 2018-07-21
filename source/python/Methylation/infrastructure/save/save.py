import numpy as np


def save_enet_top(fn, names, vals):
    info = np.zeros(len(names), dtype=[('var1', 'U50'), ('var2', int)])
    fmt = "%s %d"
    info['var1'] = list(names)
    info['var2'] = list(vals)
    np.savetxt(fn, info, fmt=fmt)

def save_names(fn, names):
    info = np.zeros(len(names), dtype=[('var1', 'U50')])
    fmt = "%s"
    info['var1'] = list(names)
    np.savetxt(fn, info, fmt=fmt)

def save_R2(fn, counts, R2s):
    info = np.zeros(len(counts), dtype=[('var1', int), ('var2', float)])
    fmt = "%d %0.18e"
    info['var1'] = counts
    info['var2'] = R2s
    np.savetxt(fn, info, fmt=fmt)

def save_linreg_top(fn, genes, p_vals, r_vals, slopes, intercepts):
    info = np.zeros(len(genes), dtype=[('var1', 'U50'), ('var2', float), ('var3', float), ('var4', float), ('var5', float)])
    fmt = "%s %0.18e %0.18e %0.18e %0.18e"
    info['var1'] = genes
    info['var2'] = p_vals
    info['var3'] = r_vals
    info['var4'] = slopes
    info['var5'] = intercepts
    np.savetxt(fn, info, fmt=fmt)

