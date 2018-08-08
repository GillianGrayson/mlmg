import numpy as np
from infrastructure.path import get_gene_data_path
from annotations.regular import get_dict_cpg_gene, get_dict_cpg_map_info
from infrastructure.load.cpg_data import load_cpg_data, get_non_inc_cpgs
from infrastructure.load.attributes import get_attributes

def load_gene_data(config):
    indexes = config.indexes
    genes_passed = []
    vals_passed = []
    fn = 'gene_data.txt'
    path = get_gene_data_path(config, fn)
    f = open(path)
    for line in f:
        col_vals = line.split(' ')
        gene = col_vals[0]
        vals = list(map(float, col_vals[1::]))
        vals = list(np.array(vals)[indexes])
        nans = np.isnan(vals)
        if True not in nans:
            genes_passed.append(gene)
            vals_passed.append(vals)

    return genes_passed, vals_passed

def get_raw_dict(config):
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_map = get_dict_cpg_map_info(config)
    attributes = get_attributes(config)
    cpgs, vals = load_cpg_data(config)

    cpg_non_inc = get_non_inc_cpgs(config)

    gene_raw_dict = {}
    map_dict = {}
    for id in range(0, len(cpgs)):

        curr_cpg = cpgs[id]
        curr_vals = vals[id]

        if curr_cpg not in cpg_non_inc:

            genes = dict_cpg_gene.get(curr_cpg)
            map_info = dict_cpg_map.get(curr_cpg)

            if genes is not None:
                for gene in genes:
                    if gene in gene_raw_dict:
                        for list_id in range(0, len(attributes)):
                            gene_raw_dict[gene][list_id].append(curr_vals[list_id])
                        map_dict[gene].append(int(map_info))
                    else:
                        gene_raw_dict[gene] = []
                        for list_id in range(0, len(attributes)):
                            gene_raw_dict[gene].append([curr_vals[list_id]])
                        map_dict[gene] = []
                        map_dict[gene].append(int(map_info))

    for gene in gene_raw_dict:
        raw = gene_raw_dict[gene]
        map_info = map_dict[gene]
        order = np.argsort(map_info)
        gene_raw_dict[gene] = []
        for record in raw:
            sorted_record = list(np.array(record)[order])
            gene_raw_dict[gene].append(sorted_record)

    return gene_raw_dict