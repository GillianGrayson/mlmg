from method.enet.routines import *
from attributes.categorical import get_attributes_dict
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path import get_result_path
from infrastructure.save.features import save_features
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
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
    genes_sorted = list(np.array(gene_names)[order])
    pvals_sorted = list(np.array(pvals)[order])


    metrics_sorted_np = np.asarray(list(map(np.log10, pvals_sorted))).reshape(-1, 1)
    bandwidth = estimate_bandwidth(metrics_sorted_np)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(metrics_sorted_np)
    labels_mean_shift = list(ms.labels_)
    clusters_mean_shift = clustering_order(labels_mean_shift)
    af = AffinityPropagation().fit(metrics_sorted_np)
    labels_affinity_propagation = list(af.labels_)
    clusters_affinity_prop = clustering_order(labels_affinity_propagation)

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_sorted, clusters_mean_shift, clusters_affinity_prop, pvals_sorted])


