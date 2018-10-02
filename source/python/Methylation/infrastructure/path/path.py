from config.method import Method
from config.types import *
from infrastructure.path.bop import *
from infrastructure.path.gene import *
from infrastructure.path.cpg import *
from infrastructure.path.experiment import *

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
        path = get_bop_path(config) + \
               '/' + InfoType.data.value + \
               '/' + file_name
    elif config.data_type is DataType.gene:
        path = get_gene_path(config)+ \
               '/' + InfoType.data.value + \
               '/' + file_name
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config)+ \
               '/' + InfoType.data.value + \
               '/' + file_name
    else:
        path = ''
    return path


def get_result_path(config, file_name):
    if config.data_type is DataType.bop:
        path = get_bop_path(config) + \
               '/' + InfoType.result.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    elif config.data_type is DataType.gene:
        path = get_gene_path(config) + \
               '/' + InfoType.result.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config) + \
               '/' + InfoType.result.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    else:
        path = ''
    return path


def get_param_path(config, file_name):
    if config.data_type is DataType.bop:
        path = get_bop_path(config) + \
               '/' + InfoType.param.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    elif config.data_type is DataType.gene:
        path = get_gene_path(config) + \
               '/' + InfoType.param.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    elif config.data_type is DataType.cpg:
        path = get_cpg_path(config) + \
               '/' + InfoType.param.value + \
               '/' + get_experiment_path(config) + \
               '/' + file_name
    else:
        path = ''
    return path
