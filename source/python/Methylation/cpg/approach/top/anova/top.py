import numpy as np
from config.config import *
from attributes.categorical import get_attributes_dict
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from scipy import stats


def save_top_anova(config, num_top=500):
    attributes_dict = get_attributes_dict(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpgs = list(dict_cpg_data.keys())
    vals = list(dict_cpg_data.values())

    pvals = []
    for id in range(0, len(cpgs)):
        curr_vals = vals[id]

        curr_beta_dict = {}
        for key_age in attributes_dict:
            curr_beta_dict[key_age] = list(np.asarray(curr_vals)[attributes_dict[key_age]])

        anova_res = stats.f_oneway(*curr_beta_dict.values())
        pvals.append(anova_res.pvalue)

    order = np.argsort(pvals)
    cpgs_sorted = list(np.array(cpgs)[order])
    pvals_sorted = list(np.array(pvals)[order])
    genes_sorted = []
    pvals_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        pval = pvals_sorted[id]
        genes = dict_cpg_gene.get(cpg)
        for gene in genes:
            if gene not in genes_sorted:
                genes_sorted.append(gene)
                pvals_genes.append(pval)

    cpgs_sorted = cpgs_sorted[0:num_top]
    pvals_sorted = pvals_sorted[0:num_top]

    genes_sorted = genes_sorted[0:num_top]
    pvals_genes = pvals_genes[0:num_top]

    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [cpgs_sorted, pvals_sorted])

    config.approach_gd = GeneDataType.from_cpg
    config.dt = DataType.gene
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [genes_sorted, pvals_genes])
    config.dt = DataType.cpg
