from method.linreg_mult.routines import *
from infrastructure.load.attributes import get_main_attributes
from infrastructure.load.top import load_top_gene_data
from infrastructure.file_system import *
from infrastructure.save.features import save_features

def save_simple_linreg(config, num_top=100):
    attributes = get_main_attributes(config)
    config.scenario = Scenario.approach
    gene_names, gene_vals = load_top_gene_data(config, num_top)
    config.scenario = Scenario.validation

    p_values = []
    r_values = []
    slopes = []
    intercepts = []
    for id in range(0, len(gene_names)):
        vals = gene_vals[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(vals, attributes)
        r_values.append(r_value)
        p_values.append(p_value)
        slopes.append(slope)
        intercepts.append(intercept)

    order = np.argsort(list(map(abs, r_values)))[::-1]
    p_values = list(np.array(p_values)[order])
    r_values = list(np.array(r_values)[order])
    slopes = list(np.array(slopes)[order])
    intercepts = list(np.array(intercepts)[order])
    gene_names = list(np.array(gene_names)[order])

    fn = 'metrics.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [gene_names, p_values, r_values, slopes, intercepts])

