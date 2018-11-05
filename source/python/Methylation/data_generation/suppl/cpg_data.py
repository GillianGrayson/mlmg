import numpy as np
from config.config import *
from annotations.gene import *
from infrastructure.load.cpg_data import load_cpg_data_raw
from infrastructure.path import *


def save_cpg_data(config, fn):
    cpgs, vals = load_cpg_data_raw(config)

    f = open(fn + '.txt')
    target_cpgs = f.read().splitlines()

    str_list = []
    for cpg in target_cpgs:
        curr_vals = vals[cpgs.index(cpg)]
        curr_str = cpg + ' '
        for id in range(0, len(curr_vals)):
            curr_str += curr_vals[id] + ' '
        str_list.append(curr_str)

    np.savetxt(fn + '_data.txt', str_list, fmt="%s")



config = Config(
    data_base=DataBase.GSE63347,
    geo_type=GeoType.any
)

fn = 'hannum_cpgs'
fn = get_suppl_path(config, fn)
save_cpg_data(config, fn)
