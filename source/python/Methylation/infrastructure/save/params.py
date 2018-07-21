import numpy as np

def save_params(fn, names, vals):
    info = np.zeros(len(names), dtype=[('var1', 'U50'), ('var2', 'float')])
    fmt = "%s %0.18e"
    info['var1'] = names
    info['var2'] = vals
    np.savetxt(fn, info, fmt=fmt)