from config.types import *

def get_root(config):
    root = config.fs.value
    return root

def get_path(config, file_name):
    root = config.fs.value
    db = config.db.value
    path = root + '/' + db + '/' + file_name
    return path

def get_suppl_path(config, file_name):
    root = config.fs.value
    db = config.db.value
    path = root + '/' + db + '/suppl_data/' + config.geo_type.value +  '/' + file_name
    return path

def get_result_path(config, file_name):
    root = config.fs.value

    path = root + '/' + config.db.value + '/' + 'result' + '/' + config.dt.value

    if config.dt is DataType.gene:
        if config.scenario is Scenario.validation:
            path += '/' + config.scenario.value + '/' + config.validation.value + '/' + config.validation_method.value
            if config.validation is Validation.simple:
                if config.validation_method is Method.linreg_mult:
                    path += '/' + config.geo_type.value + '/' + \
                            'approach(' + config.approach.value + ')_' + \
                            'method(' + config.approach_method.value + ')_' + \
                            'gender(' + config.gt.value + ')_' + \
                            'order(' + config.approach_gd.value + ')_' + \
                            'vals(' + config.validation_gd.value + ')'
                elif config.validation_method is Method.match:
                    path += '/' + config.gt.value + '/' + config.validation_gd.value + '/' + config.geo_type.value
                elif config.validation_method is Method.gender_vs:
                    path += '/' + config.approach_method.value + '/' + config.approach_gd.value
                    if config.approach_gd is GeneDataType.from_cpg:
                        path += ''
                    else:
                        path += '/' + config.geo_type.value
        elif config.scenario is Scenario.approach:
            path += '/' + config.scenario.value + '/' + config.approach.value + '/' + config.approach_method.value
            path += '/' + config.gt.value + '/' + config.approach_gd.value
            if config.approach is Approach.top:
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
            elif config.approach is Approach.clustering:
                if config.approach_gd is GeneDataType.from_cpg:
                    print('TODO')
                else:
                    print('TODO')
    elif config.dt is DataType.cpg:
        if config.scenario is Scenario.validation:
            path += '/' + config.scenario.value + '/' + config.validation.value + '/' + config.validation_method.value + '/' + config.gt.value
            if config.validation is Validation.simple:
                path += ''
        elif config.scenario is Scenario.approach:
            path += '/' + config.scenario.value + '/' + config.approach.value + '/' + config.approach_method.value + '/' + config.gt.value
            if config.approach_method is Approach.top:
                path += ''
            elif config.approach_method is Approach.clustering:
                path += ''

    path += '/' + file_name
    return path

def get_gene_data_path(config, file_name):
    root = config.fs.value
    path = root + '/' + config.db.value + '/' + 'gene_data'

    if config.scenario is Scenario.validation:
        path += '/' + config.validation_gd.value
    elif config.scenario is Scenario.approach:
        path += '/' + config.approach_gd.value

    path += '/' + config.geo_type.value
    path += '/' + file_name
    return path

def get_param_path(config, file_name):
    root = config.fs.value
    path = root + '/' + config.db.value + '/' + 'param' + '/' + config.dt.value

    if config.dt is DataType.gene:
        if config.scenario is Scenario.validation:
            path += '/' + config.scenario.value + '/' + config.validation.value + '/' + config.validation_method.value
            if config.validation is Validation.simple:
                path += ''
        elif config.scenario is Scenario.approach:
            path += '/' + config.scenario.value + '/' + config.approach.value + '/' + config.approach_method.value
            path += '/' + config.gt.value
            if config.approach is Approach.top:
                path += '/' + config.approach_gd.value
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
            elif config.approach is Approach.clustering:
                path += '/' + config.approach_gd.value
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
    elif config.dt is DataType.cpg:
        if config.scenario is Scenario.validation:
            path += '/' + config.scenario.value + '/' + config.validation.value + '/' + config.validation_method.value + '/' + config.gt.value
            if config.validation is Validation.simple:
                path += ''
        elif config.scenario is Scenario.approach:
            path += '/' + config.scenario.value + '/' + config.approach.value + '/' + config.approach_method.value + '/' + config.gt.value
            if config.approach is Approach.top:
                path += ''
            elif config.approach is Approach.clustering:
                path += ''

    path += '/' + file_name
    return path