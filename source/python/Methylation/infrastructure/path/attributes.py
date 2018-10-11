def get_attributes_path(config):
    path = config.disease.value + \
           '/' + config.gender.value
    return path
