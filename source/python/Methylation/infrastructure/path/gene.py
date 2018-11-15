def get_gene_path(config):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'gene_data' + \
           '/' + config.cross_reactive.value + \
           '/' + config.snp.value + \
           '/' + config.chromosome_type.value + \
           '/' + config.geo_type.value + \
           '/' + config.gene_data_type.value
    return path