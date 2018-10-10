def get_solution_path(config):
    path = config.scenario.value + \
           '/' + config.approach.value + \
           '/' + config.method.value
    return path
