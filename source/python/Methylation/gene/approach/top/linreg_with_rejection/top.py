from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
from scipy import stats


def save_top_linreg_with_rejection(config, part=0.05, is_clustering=False):
    attributes = get_attributes(config)
    genes, vals = load_gene_data(config)
    part_int = int(part * len(attributes))

    p_values = []
    r_values = []
    slopes = []
    intercepts = []
    for id in range(0, len(genes)):
        val = vals[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(attributes, val)

        diffs = []
        for p_id in range(0, len(attributes)):
            curr_x = attributes[p_id]
            curr_y = val[p_id]
            pred_y = slope * curr_x + intercept
            diffs.append(abs(pred_y - curr_y))

        order = np.argsort(list(map(abs, diffs)))[::-1]
        bad_ids = order[0:part_int]
        good_ids = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int)
        good_ids = list(set(good_ids) - set(bad_ids))
        good_ids.sort()

        val_good = list(np.array(val)[good_ids])
        attributes_good = list(np.array(attributes)[good_ids])

        slope, intercept, r_value, p_value, std_err = stats.linregress(attributes_good, val_good)

        r_values.append(r_value)
        p_values.append(p_value)
        slopes.append(slope)
        intercepts.append(intercept)

    order_mean = np.argsort(list(map(abs, r_values)))[::-1]
    p_values_sorted = list(np.array(p_values)[order_mean])
    r_values_sorted = list(np.array(r_values)[order_mean])
    slopes_sorted = list(np.array(slopes)[order_mean])
    intercepts_sorted = list(np.array(intercepts)[order_mean])
    genes_sorted = list(np.array(genes)[order_mean])

    features = [
        genes_sorted,
        r_values_sorted,
        p_values_sorted,
        slopes_sorted,
        intercepts_sorted
    ]
    if is_clustering:
        metrics_sorted_np = np.asarray(list(map(abs, r_values_sorted))).reshape(-1, 1)
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
