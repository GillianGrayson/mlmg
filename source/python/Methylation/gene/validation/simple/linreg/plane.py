from method.linreg_mult.routines import *
from config.types import *
from infrastructure.load.attributes import get_main_attributes
from infrastructure.load.top import load_top_gene_data, load_top_gene_names, load_top_gene_vals
from infrastructure.file_system import *
from infrastructure.save.features import save_features

def save_plane_linreg(config, num_top=100, gd_type_x=GeneDataType.mean, gd_type_y=GeneDataType.mean):
    attributes = get_main_attributes(config)
    config.scenario = Scenario.approach
    gene_names = load_top_gene_names(config, num_top)
    config.scenario = Scenario.validation
    config.approach_gd = gd_type_x
    gene_vals_x = load_top_gene_vals(config, gene_names)
    config.approach_gd = gd_type_y
    gene_vals_y = load_top_gene_vals(config, gene_names)

    p_values_x = []
    r_values_x = []
    p_values_y = []
    r_values_y = []
    for id in range(0, len(gene_names)):
        vals_main = gene_vals_x[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(vals_main, attributes)
        p_values_x.append(p_value)
        r_values_x.append(r_value)

        vals_aux = gene_vals_y[id]
        slope, intercept, r_value, p_value, std_err = stats.linregress(vals_aux, attributes)
        p_values_y.append(p_value)
        r_values_y.append(r_value)

    fn = 'plane.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [gene_names, r_values_x, r_values_y])

