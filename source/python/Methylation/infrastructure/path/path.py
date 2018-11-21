from infrastructure.path.bop import *
from infrastructure.path.gene import *
from infrastructure.path.cpg import *
from infrastructure.path.experiment import *
from infrastructure.path.attributes import *
from config.types.common import *
import os.path


def get_root(config):
    return config.fs.value


def get_origin_path(config, file_name):
    path = config.fs.value + \
           '/'  + file_name
    return path


def get_path(config, file_name):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + file_name
    return path


def get_suppl_path(config, file_name):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'suppl_data' + \
           '/' + config.geo_type.value + \
           '/' + file_name
    return path


def get_data_path(config, file_name):
    if config.data_type is DataType.bop:
        path = get_bop_path(config)
    elif config.data_type is DataType.gene:
        path = get_gene_path(config)
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config)
    else:
        path = ''

    path += '/' + InfoType.data.value

    if not os.path.exists(path):
        os.makedirs(path)

    path += '/' + file_name

    return path


def get_result_path(config, file_name):
    if config.data_type is DataType.bop:
        path = get_bop_path(config)
    elif config.data_type is DataType.gene:
        path = get_gene_path(config)
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config)
    else:
        path = ''

    path += '/' + InfoType.result.value + \
            '/' + get_experiment_path(config) + \
            '/' + get_attributes_path(config)

    if not os.path.exists(path):
        os.makedirs(path)

    path += '/' + file_name

    return path


def get_param_path(config, file_name):
    if config.data_type is DataType.bop:
        path = get_bop_path(config)
    elif config.data_type is DataType.gene:
        path = get_gene_path(config)
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config)
    else:
        path = ''

    path += '/' + InfoType.result.value + \
            '/' + get_experiment_path(config) + \
            '/' + get_attributes_path(config)

    if not os.path.exists(path):
        os.makedirs(path)

    path += '/' + file_name

    return path
