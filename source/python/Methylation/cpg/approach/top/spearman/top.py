import numpy as np
from infrastructure.load.attributes import get_attributes
from infrastructure.load.cpg_data import load_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from config.types import *
from scipy import stats


def save_top_spearman(config, num_top=500):
    attributes = get_attributes(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    cpgs, vals = load_cpg_data(config)

    rhos = []
    for id in range(0, len(cpgs)):
        curr_vals = vals[id]
        rho, pval = stats.spearmanr(attributes, curr_vals)
        rhos.append(rho)

    order = np.argsort(list(map(abs, rhos)))[::-1]
    cpgs_sorted = list(np.array(cpgs)[order])
    rhos_sorted = list(np.array(rhos)[order])

    genes_sorted = []
    rhos_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        rho = rhos_sorted[id]
        genes = dict_cpg_gene.get(cpg)
        for gene in genes:
            genes_sorted.append(gene)
            rhos_genes.append(rho)

    cpgs_sorted = cpgs_sorted[0:num_top]
    rhos_sorted = rhos_sorted[0:num_top]

    genes_sorted = genes_sorted[0:num_top]
    rhos_genes = rhos_genes[0:num_top]

    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [cpgs_sorted, rhos_sorted])

    config.approach_gd = GeneDataType.from_cpg
    config.dt = DataType.gene
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [genes_sorted, rhos_genes])
    config.dt = DataType.cpg
