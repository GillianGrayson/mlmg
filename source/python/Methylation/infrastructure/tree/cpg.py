import os
from config.types import DataType
from infrastructure.directories.common import get_chromosome_types
from infrastructure.directories.cpg import get_dna_regions
from infrastructure.directories.experiment import get_diseases, get_genders, get_scenarios, get_approaches, get_methods
from infrastructure.directories.info_type import get_info_types
from infrastructure.path import get_path
from infrastructure.tree.routines import create_fs


def create_cpg_data_tree(config):
    config.dt = DataType.cpg

    path = get_path(config, '') + 'cpg_data'
    if not os.path.exists(path):
        os.makedirs(path)

    # Common levels
    chromosome_types = get_chromosome_types(config)

    # CPG levels
    dna_regions = get_dna_regions(config)

    # Info type levels
    info_types = get_info_types(config)

    # Experiment levels
    diseases = get_diseases(config)
    genders = get_genders(config)
    scenarios = get_scenarios(config)
    approaches = get_approaches(config)
    methods = get_methods(config)

    create_fs([chromosome_types,
               dna_regions,
               info_types,
               diseases,
               genders,
               scenarios,
               approaches,
               methods],
              ['chromosome_types',
               'dna_regions',
               'info_types',
               'diseases',
               'genders',
               'scenarios',
               'approaches',
               'methods'],
              path)