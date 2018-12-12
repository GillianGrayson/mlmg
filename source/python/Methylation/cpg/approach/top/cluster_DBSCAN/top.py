import statsmodels.api as sm
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from sklearn.cluster import DBSCAN
from sklearn import metrics
from annotations.cpg import *
from annotations.gene import get_dict_cpg_gene
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.save.features import save_features
from method.clustering.order import *

import pandas as pd


def save_top_cluster_DBSCAN(config):
    attributes = get_attributes(config)
    dict_cpg_data = load_dict_cpg_data(config)
    approved_cpgs = get_approved_cpgs(config)

    attributes_normalized = [float(x) / float(max(attributes)) for x in attributes]

    if config.method_params is None or not bool(config.method_params):
        config.method_params = {}
        config.method_params['eps'] = 0.1
        config.method_params['min_samples'] = int(max(1, np.floor(len(attributes) * 0.01)))

    cpg_names_passed = []
    estimated_number_of_clusters = []
    estimated_of_noise_points = []

    num_passed = 0

    for cpg in approved_cpgs:

        if cpg in dict_cpg_data:

            values = dict_cpg_data[cpg]

            X = np.array([attributes_normalized, values]).T
            db = DBSCAN(eps=config.method_params['eps'], min_samples=config.method_params['min_samples']).fit(X)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_
            number_of_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            noise_points = list(labels).count(-1)

            cpg_names_passed.append(cpg)
            estimated_number_of_clusters.append(number_of_clusters)
            estimated_of_noise_points.append(noise_points)

            num_passed += 1

        if num_passed % config.print_rate == 0:
            print('cpg_id: ' + str(num_passed))

    order = np.argsort(list(map(abs, estimated_number_of_clusters)))[::-1]
    cpgs_sorted = list(np.array(cpg_names_passed)[order])
    estimated_number_of_clusters_sorted = list(np.array(estimated_number_of_clusters)[order])
    estimated_of_noise_points_sorted = list(np.array(estimated_of_noise_points)[order])

    features = [
        cpgs_sorted,
        estimated_number_of_clusters_sorted,
        estimated_of_noise_points_sorted,
    ]

    suffix = get_method_suffix(config.method_params)

    fn = 'top' + suffix + '.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)

    features_dict = dict()
    features_dict_keys = ['name'] + get_method_metrics(config.method)
    for feature_id in range(0, len(features_dict_keys)):
        features_dict[features_dict_keys[feature_id]] = features[feature_id]
    df = pd.DataFrame(features_dict)
    file_xls = 'top' + suffix + '.xlsx'
    file_xls = get_result_path(config, file_xls)
    writer = pd.ExcelWriter(file_xls, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()
