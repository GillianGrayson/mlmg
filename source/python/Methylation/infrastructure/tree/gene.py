import os

from config.types import DataType
from infrastructure.directories.common import get_chromosome_types
from infrastructure.directories.experiment import get_diseases, get_genders, get_scenarios, get_approaches, get_methods
from infrastructure.directories.gene import get_geo_types, get_gene_data_types
from infrastructure.directories.info_type import get_info_types
from infrastructure.path import get_path
from infrastructure.tree.routines import create_fs


def create_gene_data_tree(config):
    config.dt = DataType.gene

    path = get_path(config, '') + 'gene_data'
    if not os.path.exists(path):
        os.makedirs(path)

    # Common levels
    chromosome_types = get_chromosome_types(config)

    # GENE levels
    geo_types = get_geo_types(config)
    gene_data_types = get_gene_data_types(config)

    # Info type levels
    info_types = get_info_types(config)

    # Experiment levels
    diseases = get_diseases(config)
    genders = get_genders(config)
    scenarios = get_scenarios(config)
    approaches = get_approaches(config)
    methods = get_methods(config)

    create_fs([chromosome_types,
               geo_types,
               gene_data_types,
               info_types,
               diseases,
               genders,
               scenarios,
               approaches,
               methods],
              ['chromosome_types',
               'geo_types',
               'gene_data_types',
               'info_types',
               'diseases',
               'genders',
               'scenarios',
               'approaches',
               'methods'],
              path)