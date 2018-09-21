from infrastructure.path import *


def load_bop_cpg_dict(config):

    dict_bop_cpg = {}

    for cpg_class in config.cpg_classes:
        fn = 'dict_bop_cpg.txt'
        fn = get_bop_data_path(config, cpg_class, fn)

        f = open(fn)

        for line in f:
            line_list = line.split(' ')
            bop = line_list[0]
            cpgs = list(map(str, line_list[1::]))
            cpgs = list(map(str.strip, cpgs))
            if bop in dict_bop_cpg:
                dict_bop_cpg[bop].extend(cpgs)
            else:
                dict_bop_cpg[bop] = cpgs

        f.close()

    return dict_bop_cpg
