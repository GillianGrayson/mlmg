from config.config import *
from annotations.gene import *
from infrastructure.load.cpg_data import load_dict_cpg_data


def save_cpg_by_gene(config, fn):
    gene_cpg_dict = get_dict_gene_cpg(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpgs = list(dict_cpg_data.keys())
    vals = list(dict_cpg_data.values())

    f = open(fn + '.txt')
    target_genes = f.read().splitlines()

    str_list = []
    for gene in target_genes:
        gene_cpgs = gene_cpg_dict[gene]
        for gene_cpg in gene_cpgs:
            if gene_cpg in cpgs:
                curr_vals = vals[cpgs.index(gene_cpg)]
                curr_str = gene + ' ' + gene_cpg
                for id in range(0, len(curr_vals)):
                    curr_str += (' ' + str(format(curr_vals[id], '0.8e')))
                str_list.append(curr_str)

    np.savetxt(fn + '_cpgs.txt', str_list, fmt="%s")



config = Config(
    db=DataBase.GSE52588,
    geo=GeoType.islands_shores
)

fn = 'top_thr_corr_graph_not_in_list'
fn = get_suppl_path(config, fn)
save_cpg_by_gene(config, fn)
