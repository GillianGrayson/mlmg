from method.enet.routines import *
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *


def save_top_moment(config):
    gene_names, gene_values = load_gene_data(config)

    means = []
    stds = []
    for id in range(0, len(gene_names)):
        values = gene_values[id]
        means.append(np.mean(values))
        stds.append(np.std(values))

    order = np.argsort(list(map(abs, means)))[::-1]
    genes_sorted = list(np.array(gene_names)[order])
    means_sorted = list(np.array(means)[order])
    stds_sorted = list(np.array(stds)[order])

    features = [
        genes_sorted,
        means_sorted,
        stds_sorted
    ]
    if config.is_clustering:
        metrics_sorted_np = np.asarray(list(map(abs, means_sorted))).reshape(-1, 1)
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
