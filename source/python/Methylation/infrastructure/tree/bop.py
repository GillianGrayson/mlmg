import os
from config.types import *
from infrastructure.directories.bop import get_class_types
from infrastructure.directories.common import get_chromosome_types
from infrastructure.directories.experiment import get_diseases, get_genders, get_scenarios, get_approaches, get_methods
from infrastructure.directories.info_type import get_info_types
from infrastructure.path.path import get_path
from infrastructure.tree.routines import create_fs


def create_bop_data_tree(config):
    config.dt = DataType.bop

    path = get_path(config, '') + 'bop_data'
    if not os.path.exists(path):
        os.makedirs(path)

    # Common levels
    chromosome_types = get_chromosome_types(config)

    # BOP levels
    class_types = get_class_types(config)

    # Info type levels
    info_types = [InfoType.result.value, InfoType.param.value]
    info_type_data = [InfoType.data.value]

    # Experiment levels
    diseases = get_diseases(config)
    genders = get_genders(config)
    scenarios = get_scenarios(config)
    approaches = get_approaches(config)
    methods = get_methods(config)

    create_fs([chromosome_types,
               class_types,
               info_type_data],
              ['chromosome_types',
               'class_types',
               'info_types'],
              path, False)

    create_fs([chromosome_types,
               class_types,
               info_types,
               diseases,
               genders,
               scenarios,
               approaches,
               methods],
              ['chromosome_types',
               'class_types',
               'info_types',
               'diseases',
               'genders',
               'scenarios',
               'approaches',
               'methods'],
              path)