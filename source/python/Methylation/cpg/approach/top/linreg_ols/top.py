import statsmodels.api as sm
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation

from annotations.cpg import *
from annotations.gene import get_dict_cpg_gene
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.save.features import save_features
from method.clustering.order import *

import pandas as pd


def save_top_linreg_ols(config):
    attributes = get_attributes(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    approved_cpgs = get_approved_cpgs(config)

    print('len(dict_cpg_data): ' + str(len(dict_cpg_data)))
    print('len(approved_cpgs): ' + str(len(approved_cpgs)))

    if not bool(config.method_params):
        config.method_params = {}
        config.method_params['outliers_limit'] = 0.0
        config.method_params['outliers_sigma'] = 0.0

    cpg_names_passed = []
    R2s = []
    intercepts = []
    slopes = []
    intercepts_std_errors = []
    slopes_std_errors = []
    intercepts_p_values = []
    slopes_p_values = []

    num_passed = 0

    for cpg in approved_cpgs:

        if cpg in dict_cpg_data:

            values = dict_cpg_data[cpg]

            x = sm.add_constant(attributes)
            results = sm.OLS(values, x).fit()

            if np.isclose(config.method_params['outliers_limit'], 0.0):
                cpg_names_passed.append(cpg)
                R2s.append(results.rsquared)
                intercepts.append(results.params[0])
                slopes.append(results.params[1])
                intercepts_std_errors.append(results.bse[0])
                slopes_std_errors.append(results.bse[1])
                intercepts_p_values.append(results.pvalues[0])
                slopes_p_values.append(results.pvalues[1])
                num_passed += 1
            else:
                slope_plus = results.params[1] + config.method_params['outliers_sigma'] * results.bse[1]
                intercept_plus = results.params[0] + config.method_params['outliers_sigma'] * results.bse[0]

                slope = results.params[1]
                intercept = results.params[0]

                max_diff = ((slope_plus * max(attributes) + intercept_plus) - (slope * max(attributes) + intercept))

                passed_ids = []
                for p_id in range(0, len(attributes)):
                    curr_x = attributes[p_id]
                    curr_y = values[p_id]
                    pred_y = results.params[1] * curr_x + results.params[0]
                    if abs(pred_y - curr_y) < max_diff:
                        passed_ids.append(p_id)

                if len(passed_ids) > np.floor(len(values) * config.method_params['outliers_limit']):

                    values_good = list(np.array(values)[passed_ids])
                    attributes_good = list(np.array(attributes)[passed_ids])

                    x = sm.add_constant(attributes_good)
                    results = sm.OLS(values_good, x).fit()

                    cpg_names_passed.append(cpg)
                    R2s.append(results.rsquared)
                    intercepts.append(results.params[0])
                    slopes.append(results.params[1])
                    intercepts_std_errors.append(results.bse[0])
                    slopes_std_errors.append(results.bse[1])
                    intercepts_p_values.append(results.pvalues[0])
                    slopes_p_values.append(results.pvalues[1])
                    num_passed += 1

        if num_passed % config.print_rate == 0:
            print('cpg_id: ' + str(num_passed))

    order = np.argsort(list(map(abs, R2s)))[::-1]
    cpgs_sorted = list(np.array(cpg_names_passed)[order])
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

    suffix = get_method_suffix(config.method_params)

    fn = 'top' + suffix + '.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)

    features_dict = dict()
    features_dict_keys = ['name'] + get_method_metrics(config.method)
    for feature_id in range(0, len(features_dict_keys)):
        features_dict[features_dict_keys[feature_id]] = features[feature_id]
    if config.is_clustering:
        features_dict['mean_shift'] = clusters_mean_shift
        features_dict['affinity_prop'] = clusters_affinity_prop
    df = pd.DataFrame(features_dict)
    file_xls = 'top' + suffix + '.xlsx'
    file_xls = get_result_path(config, file_xls)
    writer = pd.ExcelWriter(file_xls, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

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
        features = features + [
            clusters_mean_shift_genes,
            clusters_affinity_prop_genes,
        ]

    fn = 'top' + suffix + '.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)

    features_dict = dict()
    features_dict_keys = ['name'] + get_method_metrics(config.method)
    for feature_id in range(0, len(features_dict_keys)):
        features_dict[features_dict_keys[feature_id]] = features[feature_id]
    if config.is_clustering:
        features_dict['mean_shift'] = clusters_mean_shift_genes
        features_dict['affinity_prop'] = clusters_affinity_prop_genes
    df = pd.DataFrame(features_dict)
    file_xls = 'top' + suffix + '.xlsx'
    file_xls = get_result_path(config, file_xls)
    writer = pd.ExcelWriter(file_xls, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

    config.gene_data_type = gene_data_type
    config.geo_type = geo_type
    config.data_type = DataType.cpg