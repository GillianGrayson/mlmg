import numpy as np
from config.types.experiments.method import Method
from infrastructure.path.path import get_param_path


def load_params_dict(config):
    fn = 'params.txt'
    path = get_param_path(config, fn)
    params = np.loadtxt(path, dtype='U50')

    params_dict = {}
    if config.method is Method.enet:
        alpha = float(params[0][1])
        l1_ratio = float(params[1][1])
        params_dict['alpha'] = alpha
        params_dict['l1_ratio'] = l1_ratio

    return params_dict
