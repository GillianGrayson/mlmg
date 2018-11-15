import numpy as np
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from config.config import *
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
from annotations.cpg import *


def save_top_moment(config):
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpg_names = list(dict_cpg_data.keys())
    cpg_values = list(dict_cpg_data.values())
    approved_cpgs = get_approved_cpgs(config)

    cpg_names_passed = []
    means = []
    stds = []
    for id in range(0, len(cpg_names)):
        cpg = cpg_names[id]
        if cpg in approved_cpgs:
            cpg_names_passed.append(cpg)
            values = cpg_values[id]
            means.append(np.mean(values))
            stds.append(np.std(values))

            if id % config.print_rate == 0:
                print('cpg_id: ' + str(id))

    order = np.argsort(list(map(abs, means)))[::-1]
    cpgs_sorted = list(np.array(cpg_names_passed)[order])
    means_sorted = list(np.array(means)[order])
    stds_sorted = list(np.array(stds)[order])

    features = [
        cpgs_sorted,
        means_sorted,
        stds_sorted
    ]

    clusters_mean_shift = []
    clusters_affinity_prop = []
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
        features =  features + [
            clusters_mean_shift,
            clusters_affinity_prop,
        ]
        fn = 'top_with_clustering.txt'
    else:
        fn = 'top.txt'

    fn = get_result_path(config, fn)
    save_features(fn, features)

    genes_sorted = []
    cpgs_genes = []
    means_genes = []
    stds_genes = []
    clusters_mean_shift_genes = []
    clusters_affinity_prop_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        mean = means_sorted[id]
        std = stds_sorted[id]
        if config.is_clustering:
            cluster_mean_shift = clusters_mean_shift[id]
            cluster_affinity_prop = clusters_affinity_prop[id]
        else:
            cluster_mean_shift = []
            cluster_affinity_prop = []
        if cpg in dict_cpg_gene:
            genes = dict_cpg_gene.get(cpg)
            for gene in genes:
                if gene not in genes_sorted:
                    genes_sorted.append(gene)
                    cpgs_genes.append(cpg)
                    means_genes.append(mean)
                    stds_genes.append(std)
                    if config.is_clustering:
                        clusters_mean_shift_genes.append(cluster_mean_shift)
                        clusters_affinity_prop_genes.append(cluster_affinity_prop)

    config.data_type = DataType.gene
    gene_data_type = config.gene_data_type
    geo_type = config.geo_type
    config.gene_data_type = GeneDataType.from_cpg
    config.geo_type = GeoType.from_cpg

    features = [
        genes_sorted,
        cpgs_genes,
        means_genes,
        stds_genes
    ]

    if config.is_clustering:
        fn = 'top_with_clustering.txt'
        features = features + [
            clusters_mean_shift_genes,
            clusters_affinity_prop_genes,
        ]
    else:
        fn = 'top.txt'

    fn = get_result_path(config, fn)
    save_features(fn, features)

    config.gene_data_type = gene_data_type
    config.geo_type = geo_type
    config.data_type = DataType.cpg