def get_cpg_path(config):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'gene_data' + \
           '/' + config.chromosome_type.value + \
           '/' + config.dna_region.value
    return path