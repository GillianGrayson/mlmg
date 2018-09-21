import numpy as np
from config.config import *
from annotations.regular import *
from infrastructure.load.cpg_data import load_cpg_data
from infrastructure.path import *


def save_cpg_data(config, fn):
    cpgs, vals = load_cpg_data(config)

    f = open(fn + '.txt')
    target_cpgs = f.read().splitlines()

    str_list = []
    for cpg in target_cpgs:
        curr_vals = vals[cpgs.index(cpg)]
        curr_str = cpg
        for id in range(0, len(curr_vals)):
            curr_str += (' ' + str(format(curr_vals[id], '0.8e')))
        str_list.append(curr_str)

    np.savetxt(fn + '_cpgs.txt', str_list, fmt="%s")



config = Config(
    db=DataBaseType.GSE52588,
    geo=GeoType.any
)

fn = 'horvatz_cpgs'
fn = get_suppl_path(config, fn)
save_cpg_data(config, fn)
