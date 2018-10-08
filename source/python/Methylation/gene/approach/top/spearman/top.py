from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from scipy import stats
import math
from method.clustering.order import *
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation


def save_top_spearman(config, is_clustering=False):
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

    features = [
        genes_sorted,
        rhos_sorted
    ]
    if is_clustering:
        metrics_sorted_np = np.asarray(list(map(abs, rhos_sorted))).reshape(-1, 1)
        bandwidth = estimate_bandwidth(metrics_sorted_np)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(metrics_sorted_np)
        labels_mean_shift = list(ms.labels_)
        clusters_mean_shift = clustering_order(labels_mean_shift)
        af = AffinityPropagation().fit(metrics_sorted_np)
        labels_affinity_propagation = list(af.labels_)
        clusters_affinity_prop = clustering_order(labels_affinity_propagation)
        features = features + [
            clusters_mean_shift,
            clusters_affinity_prop
        ]
        fn = 'top_with_clustering.txt'
    else:
        fn = 'top.txt'

    fn = get_result_path(config, fn)
    save_features(fn, features)
