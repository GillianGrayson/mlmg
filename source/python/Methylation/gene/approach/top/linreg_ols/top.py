from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
import statsmodels.api as sm


def save_top_linreg_ols(config):
    attributes = get_attributes(config)
    gene_names, gene_values = load_gene_data(config)

    R2s = []
    intercepts = []
    slopes = []
    intercepts_std_errors = []
    slopes_std_errors = []
    intercepts_p_values = []
    slopes_p_values = []
    for id in range(0, len(gene_names)):
        values = gene_values[id]
        x = sm.add_constant(attributes)
        results = sm.OLS(values, x).fit()
        R2s.append(results.rsquared)
        intercepts.append(results.params[0])
        slopes.append(results.params[1])
        intercepts_std_errors.append(results.bse[0])
        slopes_std_errors.append(results.bse[1])
        intercepts_p_values.append(results.pvalues[0])
        slopes_p_values.append(results.pvalues[1])

    order = np.argsort(list(map(abs, R2s)))[::-1]
    genes_sorted = list(np.array(gene_names)[order])
    R2s_sorted = list(np.array(R2s)[order])
    intercepts_sorted = list(np.array(intercepts)[order])
    slopes_sorted = list(np.array(slopes)[order])
    intercepts_std_errors_sorted = list(np.array(intercepts_std_errors)[order])
    slopes_std_errors_sorted = list(np.array(slopes_std_errors)[order])
    intercepts_p_values_sorted = list(np.array(intercepts_p_values)[order])
    slopes_p_values_sorted = list(np.array(slopes_p_values)[order])

    features = [
        genes_sorted,
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
        features = features + [
            clusters_mean_shift,
            clusters_affinity_prop
        ]
        fn = 'top_with_clustering.txt'
    else:

        fn = 'top.txt'

    fn = get_result_path(config, fn)
    save_features(fn, features)
