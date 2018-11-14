def get_cpg_path(config):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'cpg_data' + \
           '/' + config.cross_reactive.value + \
           '/' + config.snp.value + \
           '/' + config.chromosome_type.value + \
           '/' + config.dna_region.value
    return path