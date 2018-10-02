def get_gene_path(config):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'gene_data' + \
           '/' + config.chromosome_type.value + \
           '/' + config.geo_type.value + \
           '/' + config.gene_data_type.value
    return path