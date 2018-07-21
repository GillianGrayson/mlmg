from method.enet.routines import *
from infrastructure.load.attributes import get_main_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.file_system import get_result_path
from infrastructure.save.features import save_features
from scipy import stats


def save_top_anova(config, shift=5):
    attributes = get_main_attributes(config)
    gene_names, gene_vals = load_gene_data(config)

    min_atr = min(attributes)
    max_atr = max(attributes)

    min_atr = int(min_atr / shift) * shift
    max_atr = (int(max_atr / shift) + 1) * shift
    age_dict = {}
    for age_id in range(0, len(attributes)):
        age = attributes[age_id]
        key = int((age - min_atr) / shift)
        if key in age_dict:
            age_dict[key].append(age_id)
        else:
            age_dict[key] = [age_id]

    pvals = []
    for id in range(0, len(gene_names)):

        vals = gene_vals[id]

        vals_dict = {}
        for key_age in age_dict:
            vals_dict[key_age] = list(np.asarray(vals)[age_dict[key_age]])

        anova_mean = stats.f_oneway(*vals_dict.values())
        pvals.append(anova_mean.pvalue)

    order = np.argsort(pvals)
    genes_opt = list(np.array(gene_names)[order])
    pvals_opt = list(np.array(pvals)[order])

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_opt, pvals_opt])

