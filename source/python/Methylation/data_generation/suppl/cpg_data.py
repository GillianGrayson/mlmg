from config.config import *
from annotations.gene import *
from infrastructure.load.cpg_data import load_cpg_data_raw


def save_cpg_data(config, fn):
    cpgs, vals = load_cpg_data_raw(config)

    f = open(fn + '.txt')
    target_cpgs = f.read().splitlines()

    str_list = []
    for cpg in target_cpgs:
        if cpg in cpgs:
            curr_vals = vals[cpgs.index(cpg)]
            curr_str = cpg + ' ' + ' '.join(curr_vals)
            str_list.append(curr_str)

    np.savetxt(fn + '_data.txt', str_list, fmt="%s")



config = Config(
    data_base=DataBase.GSE40279,
    geo_type=GeoType.any
)

fn = 'GSE40279_4classes_new_list_best_cpg_from_1000_top_genes'
fn = get_suppl_path(config, fn)
save_cpg_data(config, fn)
