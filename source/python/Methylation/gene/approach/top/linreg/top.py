from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path import get_result_path
from infrastructure.save.features import save_features
from scipy import stats


def save_top_linreg(config):
    attributes = get_attributes(config)
    genes, vals = load_gene_data(config)

    p_values = []
    r_values = []
    slopes = []
    intercepts = []
    for id in range(0, len(genes)):
        val = vals[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(val, attributes)
        r_values.append(r_value)
        p_values.append(p_value)
        slopes.append(slope)
        intercepts.append(intercept)

    order_mean = np.argsort(list(map(abs, r_values)))[::-1]
    p_values_opt = list(np.array(p_values)[order_mean])
    r_values_opt = list(np.array(r_values)[order_mean])
    slopes_opt = list(np.array(slopes)[order_mean])
    intercepts_opt = list(np.array(intercepts)[order_mean])
    genes_opt = list(np.array(genes)[order_mean])

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_opt, r_values_opt, p_values_opt, slopes_opt, intercepts_opt])
