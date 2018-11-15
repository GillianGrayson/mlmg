from config.types.annotations import GeneDataType, GeoType


def get_geo_types(config):
    geo_types = [GeoType.islands_shores.value]
    return geo_types

def get_gene_data_types(config):
    gene_data_types = [GeneDataType.mean.value]
    return gene_data_types