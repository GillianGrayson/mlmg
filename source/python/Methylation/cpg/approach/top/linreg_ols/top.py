import numpy as np
from infrastructure.load.attributes import get_attributes
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from config.config import *
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
import statsmodels.api as sm


def save_top_linreg_ols(config):
    attributes = get_attributes(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpg_names = list(dict_cpg_data.keys())
    cpg_values = list(dict_cpg_data.values())

    R2s = []
    intercepts = []
    slopes = []
    intercepts_std_errors = []
    slopes_std_errors = []
    intercepts_p_values = []
    slopes_p_values = []
    for id in range(0, len(cpg_names)):
        values = cpg_values[id]
        x = sm.add_constant(attributes)
        results = sm.OLS(values, x).fit()
        R2s.append(results.rsquared)
        intercepts.append(results.params[0])
        slopes.append(results.params[1])
        intercepts_std_errors.append(results.bse[0])
        slopes_std_errors.append(results.bse[1])
        intercepts_p_values.append(results.pvalues[0])
        slopes_p_values.append(results.pvalues[1])

        if id % config.print_rate == 0:
            print('cpg_id: ' + str(id))

    order = np.argsort(list(map(abs, R2s)))[::-1]
    cpgs_sorted = list(np.array(cpg_names)[order])
    R2s_sorted = list(np.array(R2s)[order])
    intercepts_sorted = list(np.array(intercepts)[order])
    slopes_sorted = list(np.array(slopes)[order])
    intercepts_std_errors_sorted = list(np.array(intercepts_std_errors)[order])
    slopes_std_errors_sorted = list(np.array(slopes_std_errors)[order])
    intercepts_p_values_sorted = list(np.array(intercepts_p_values)[order])
    slopes_p_values_sorted = list(np.array(slopes_p_values)[order])

    clusters_mean_shift = []
    clusters_affinity_prop = []
    features = [
        cpgs_sorted,
        R2s_sorted,
        intercepts_sorted,
        slopes_sorted,
        intercepts_std_errors_sorted,
        slopes_std_errors_sorted,
        intercepts_p_values_sorted,
        slopes_p_values_sorted
    ]
    if config.is_clustering:
        metrics_sorted_np = np.asarray(list(map(abs, R2s_sorted))).reshape(-1, 1)
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
    R2s_genes = []
    intercepts_genes = []
    slopes_genes = []
    intercepts_std_errors_genes = []
    slopes_std_errors_genes = []
    intercepts_p_values_genes = []
    slopes_p_values_genes = []
    clusters_mean_shift_genes = []
    clusters_affinity_prop_genes = []
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        R2 = R2s_sorted[id]
        intercept = intercepts_sorted[id]
        slope = slopes_sorted[id]
        intercepts_std_error = intercepts_std_errors_sorted[id]
        slopes_std_error = slopes_std_errors_sorted[id]
        intercepts_p_value = intercepts_p_values_sorted[id]
        slopes_p_value = slopes_p_values_sorted[id]
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
                    R2s_genes.append(R2)
                    intercepts_genes.append(intercept)
                    slopes_genes.append(slope)
                    intercepts_std_errors_genes.append(intercepts_std_error)
                    slopes_std_errors_genes.append(slopes_std_error)
                    intercepts_p_values_genes.append(intercepts_p_value)
                    slopes_p_values_genes.append(slopes_p_value)
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
        R2s_genes,
        intercepts_genes,
        slopes_genes,
        intercepts_std_errors_genes,
        slopes_std_errors_genes,
        intercepts_p_values_genes,
        slopes_p_values_genes,
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