from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path import *
from infrastructure.save.features import save_features


def save_params_enet(config, num_folds=10):
    attributes = get_attributes(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpgs = list(dict_cpg_data.keys())
    vals = list(dict_cpg_data.values())

    param_names, param_values = get_enet_params(attributes, vals, num_folds)

    fn = 'params.txt'
    fn = get_param_path(config, fn)
    save_features(fn, [param_names, param_values])

