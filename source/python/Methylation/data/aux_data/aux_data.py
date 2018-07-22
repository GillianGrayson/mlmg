import numpy as np
from annotation.regular import *
from infrastructure.load.cpg_data import load_cpg_data
from infrastructure.file_system import *


def save_cpg_by_gene_list(config, fn):
    gene_cpg_dict = get_dict_gene_cpg(config)
    cpgs, vals = load_cpg_data(config)

    fn = get_aux_path(config, fn)
    f = open(fn + '.txt')
    target_genes = f.read().splitlines()

    str_list = []
    for gene in target_genes:
        gene_cpgs = gene_cpg_dict[gene]
        for gene_cpg in gene_cpgs:
            if gene_cpg in cpgs:
                index = cpgs.index(gene_cpg)
                curr_vals = vals[cpgs.index(gene_cpg)]
                curr_str = gene + ' ' + gene_cpg
                for id in range(0, len(curr_vals)):
                    curr_str += (' ' + str(format(curr_vals[id], '0.8e')))
                str_list.append(curr_str)

    np.savetxt(fn + '_cpgs.txt', str_list, fmt="%s")
