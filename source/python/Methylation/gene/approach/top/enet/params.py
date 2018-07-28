from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.file_system import *
from infrastructure.save.features import save_features


def save_params_enet(config, num_folds=10):

    attributes = get_attributes(config)
    genes_passed, vals_passed = load_gene_data(config)

    param_names, param_values = get_enet_params(vals_passed, attributes, num_folds)

    fn = 'params.txt'
    fn = get_param_path(config, fn)
    save_features(fn, [param_names, param_values])






