def get_bop_path(config):
    path = config.fs.value + \
           '/' + config.data_base.value + \
           '/' + 'bop_data' + \
           '/' + config.cross_reactive.value + \
           '/' + config.snp.value + \
           '/' + config.chromosome_type.value + \
           '/' + config.class_type.value
    return path
