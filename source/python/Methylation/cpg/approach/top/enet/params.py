from method.enet.routines import *
from infrastructure.load.attributes import get_main_attributes
from infrastructure.load.cpg_data import load_cpg_data
from infrastructure.file_system import *
from infrastructure.save.features import save_features


def save_params_enet(config, num_folds=10):
    attributes = get_main_attributes(config)
    cpgs, vals = load_cpg_data(config)

    param_names, param_values = get_enet_params(attributes, vals, num_folds)

    fn = 'params.txt'
    fn = get_param_path(config, fn)
    save_features(fn, [param_names, param_values])

