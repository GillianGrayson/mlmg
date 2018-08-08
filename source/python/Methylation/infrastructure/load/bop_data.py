from infrastructure.path import *


def load_bop_cpg_dict(config):
    fn = 'dict_bop_cpg.txt'
    fn = get_bop_data_path(config, fn)

    f = open(fn)

    dict_bop_cpg = {}
    for line in f:
        line_list = line.split(' ')
        bop = line_list[0]
        cpgs = list(map(str, line_list[1::]))
        cpgs = list(map(str.strip, cpgs))
        dict_bop_cpg[bop] = cpgs

    return dict_bop_cpg
