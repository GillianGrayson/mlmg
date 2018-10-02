def get_experiment_path(config):
    path = config.disease.value + \
           '/' + config.gender.value + \
           '/' + config.scenario.value + \
           '/' + config.approach.value + \
           '/' + config.method.value
    return path