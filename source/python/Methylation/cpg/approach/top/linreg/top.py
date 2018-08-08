import numpy as np
from infrastructure.load.attributes import get_attributes
from infrastructure.load.cpg_data import load_cpg_data
from infrastructure.path import get_result_path
from infrastructure.save.features import save_features
from annotations.regular import get_dict_cpg_gene
from config.types import *
from scipy import stats


def save_top_linreg(config, num_top=500):
    attributes = get_attributes(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    cpgs, vals = load_cpg_data(config)

    slopes = []
    intercepts = []
    rvals = []
    pvals = []
    for id in range(0, len(cpgs)):
        curr_vals = vals[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(curr_vals, attributes)
        slopes.append(slope)
        intercepts.append(intercept)
        rvals.append(r_value)
        pvals.append(p_value)

    order = np.argsort(pvals)
    cpgs_sorted = list(np.array(cpgs)[order])
    pvals_sorted = list(np.array(pvals)[order])
    slopes_sorted = list(np.array(slopes)[order])
    intercepts_sorted = list(np.array(intercepts)[order])
    rvals_sorted = list(np.array(rvals)[order])

    genes_sorted = []
    pvals_genes = []
    slopes_genes = []
    intercepts_genes = []
    rvals_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        pval = pvals_sorted[id]
        slope = slopes_sorted[id]
        intercept = intercepts_sorted[id]
        rval = rvals_sorted[id]
        genes = dict_cpg_gene.get(cpg)
        for gene in genes:
            if gene not in genes_sorted:
                genes_sorted.append(gene)
                pvals_genes.append(pval)
                slopes_genes.append(slope)
                intercepts_genes.append(intercept)
                rvals_genes.append(rval)

    cpgs_sorted = cpgs_sorted[0:num_top]
    pvals_sorted = pvals_sorted[0:num_top]
    slopes_sorted = slopes_sorted[0:num_top]
    intercepts_sorted = intercepts_sorted[0:num_top]
    rvals_sorted = rvals_sorted[0:num_top]

    genes_sorted = genes_sorted[0:num_top]
    pvals_genes = pvals_genes[0:num_top]
    slopes_genes = slopes_genes[0:num_top]
    intercepts_genes = intercepts_genes[0:num_top]
    rvals_genes = rvals_genes[0:num_top]

    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [cpgs_sorted, pvals_sorted, rvals_sorted, slopes_sorted, intercepts_sorted])

    config.approach_gd = GeneDataType.from_cpg
    config.dt = DataType.gene
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [genes_sorted, pvals_genes, rvals_genes, slopes_genes, intercepts_genes])
    config.dt = DataType.cpg