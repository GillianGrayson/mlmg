from method.enet.routines import *
from method.anova.routines import get_attributes_dict
from infrastructure.load.gene_data import load_gene_data
from infrastructure.file_system import get_result_path
from infrastructure.save.features import save_features
from scipy import stats


def save_top_anova(config):
    gene_names, gene_vals = load_gene_data(config)
    attributes_dict = get_attributes_dict(config)

    pvals = []
    for id in range(0, len(gene_names)):

        vals = gene_vals[id]

        vals_dict = {}
        for key_age in attributes_dict:
            vals_dict[key_age] = list(np.asarray(vals)[attributes_dict[key_age]])

        anova_mean = stats.f_oneway(*vals_dict.values())
        pvals.append(anova_mean.pvalue)

    order = np.argsort(pvals)
    genes_opt = list(np.array(gene_names)[order])
    pvals_opt = list(np.array(pvals)[order])

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_opt, pvals_opt])


