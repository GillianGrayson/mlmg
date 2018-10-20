import numpy as np
from infrastructure.load.attributes import get_attributes
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from config.types import *
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
from scipy import stats


def save_top_linreg(config):
    attributes = get_attributes(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpgs = list(dict_cpg_data.keys())
    vals = list(dict_cpg_data.values())

    r_values = []
    p_values = []
    slopes = []
    intercepts = []
    std_errors = []
    for id in range(0, len(cpgs)):
        curr_vals = vals[id]
        slope, intercept, r_value, p_value, std_error = stats.linregress(attributes, curr_vals)
        r_values.append(r_value)
        p_values.append(p_value)
        slopes.append(slope)
        intercepts.append(intercept)
        std_errors.append(std_error)
        if id % config.print_rate == 0:
            print('cpg_id: ' + str(id))

    order = np.argsort(list(map(abs, r_values)))[::-1]
    cpgs_sorted = list(np.array(cpgs)[order])
    r_values_sorted = list(np.array(r_values)[order])
    p_values_sorted = list(np.array(p_values)[order])
    slopes_sorted = list(np.array(slopes)[order])
    intercepts_sorted = list(np.array(intercepts)[order])
    std_errors_sorted = list(np.array(std_errors)[order])

    clusters_mean_shift = []
    clusters_affinity_prop = []
    features = [
        cpgs_sorted,
        r_values_sorted,
        p_values_sorted,
        slopes_sorted,
        intercepts_sorted,
        std_errors_sorted
    ]
    if config.is_clustering:
        metrics_sorted_np = np.asarray(list(map(abs, r_values_sorted))).reshape(-1, 1)
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
    clusters_mean_shift_genes = []
    clusters_affinity_prop_genes = []
    r_values_genes = []
    p_values_genes = []
    slopes_genes = []
    intercepts_genes = []
    std_errors_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        if config.is_clustering:
            cluster_mean_shift = clusters_mean_shift[id]
            cluster_affinity_prop = clusters_affinity_prop[id]
        else:
            cluster_mean_shift = []
            cluster_affinity_prop = []
        r_value = r_values_sorted[id]
        p_value = p_values_sorted[id]
        slope = slopes_sorted[id]
        intercept = intercepts_sorted[id]
        std_error = std_errors_sorted[id]
        if cpg in dict_cpg_gene:
            genes = dict_cpg_gene.get(cpg)
            for gene in genes:
                if gene not in genes_sorted:
                    genes_sorted.append(gene)
                    cpgs_genes.append(cpg)
                    if config.is_clustering:
                        clusters_mean_shift_genes.append(cluster_mean_shift)
                        clusters_affinity_prop_genes.append(cluster_affinity_prop)
                    r_values_genes.append(r_value)
                    p_values_genes.append(p_value)
                    slopes_genes.append(slope)
                    intercepts_genes.append(intercept)
                    std_errors_genes.append(std_error)

    config.data_type = DataType.gene
    gene_data_type = config.gene_data_type
    geo_type = config.geo_type
    config.gene_data_type = GeneDataType.from_cpg
    config.geo_type = GeoType.from_cpg

    features = [
        genes_sorted,
        cpgs_genes,
        r_values_genes,
        p_values_genes,
        slopes_genes,
        intercepts_genes,
        std_errors_genes
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