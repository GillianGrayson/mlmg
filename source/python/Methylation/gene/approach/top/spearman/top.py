from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.file_system import get_result_path
from infrastructure.save.features import save_features
from scipy import stats
import math


def save_top_spearman(config):
    attributes = get_attributes(config)
    gene_names, gene_vals = load_gene_data(config)

    gene_rhos = []

    for id in range(0, len(gene_names)):

        vals = gene_vals[id]
        rhos, pvals = stats.spearmanr(attributes, vals)
        if math.isnan(rhos):
            rhos = 0.0
        gene_rhos.append(rhos)

    order = np.argsort(list(map(abs, gene_rhos)))[::-1]
    rhos_sorted = list(np.array(gene_rhos)[order])
    genes_sorted = list(np.array(gene_names)[order])

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_sorted, rhos_sorted])
