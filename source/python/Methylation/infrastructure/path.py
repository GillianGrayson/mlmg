from config.types import *

def get_root(config):
    return config.fs.value

def get_origin_path(config, file_name):
    path = config.fs.value + '/'  + file_name
    return path

def get_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/' + file_name
    return path

def get_suppl_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/suppl_data/' + config.geo_type.value +  '/' + file_name
    return path

def get_result_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/' + 'result'

    path += '/' + config.dt.value
    if config.dt is DataType.gene:

        path += '/' + config.scenario.value
        if config.scenario is Scenario.validation:

            path += '/' + config.validation.value
            if config.validation is Validation.simple:

                path += '/' + config.validation_method.value
                if config.validation_method is Method.linreg_mult:

                    path += '/' + config.approach.value + \
                            '/' + config.approach_method.value + \
                            '/' + config.gender.value + \
                            '/' + config.disease.value + \
                            '/' + config.geo_type.value + \
                            '/' + config.approach_gd.value + \
                            '/' + config.validation_gd.value

                elif config.validation_method is Method.match:

                    path += '/' + config.gender.value + \
                            '/' + config.disease.value + \
                            '/' + config.validation_gd.value + \
                            '/' + config.geo_type.value

                elif config.validation_method is Method.gender_vs:

                    path += '/' + config.approach_method.value + \
                            '/' + config.approach_gd.value

                    if config.approach_gd is GeneDataType.from_cpg:
                        path += ''
                    elif config.approach_gd is GeneDataType.from_bop:
                        path += ''
                    else:
                        path += '/' + config.geo_type.value

        elif config.scenario is Scenario.approach:

            path += '/' + config.approach.value + \
                    '/' + config.approach_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value + \
                    '/' + config.approach_gd.value

            if config.approach is Approach.top:
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                elif config.approach_gd is GeneDataType.from_bop:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
            if config.approach is Approach.bend:
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                elif config.approach_gd is GeneDataType.from_bop:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
            elif config.approach is Approach.clustering:
                path += ''

    elif config.dt is DataType.cpg:

        path += '/' + config.scenario.value
        if config.scenario is Scenario.validation:

            path += '/' + config.validation.value + \
                    '/' + config.validation_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value

            if config.validation is Validation.simple:
                path += ''

        elif config.scenario is Scenario.approach:

            if config.approach is Approach.top:

                path += '/' + config.approach.value + \
                        '/' + config.approach_method.value + \
                        '/' + config.gender.value + \
                        '/' + config.disease.value

            elif config.approach is Approach.clustering:

                path += ''

            elif config.approach is Approach.statistics:

                path += '/' + config.approach.value + \
                        '/' + config.gender.value + \
                        '/' + config.disease.value + \
                        '/' + config.geo_type.value

    elif config.dt is DataType.bop:

        path += '/' + config.scenario.value
        if config.scenario is Scenario.validation:

            path += '/' + config.validation.value + \
                    '/' + config.validation_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value

            if config.validation is Validation.simple:
                path += ''

        elif config.scenario is Scenario.approach:

            path += '/' + config.approach.value + \
                    '/' + config.approach_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value

            if config.approach_method is Approach.top:
                path += ''
            elif config.approach_method is Approach.clustering:
                path += ''

    path += '/' + file_name
    return path

def get_gene_data_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/' + 'gene_data'

    if config.scenario is Scenario.validation:
        path += '/' + config.validation_gd.value
    elif config.scenario is Scenario.approach:
        path += '/' + config.approach_gd.value

    path += '/' + config.geo_type.value
    path += '/' + file_name
    return path

def get_bop_data_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/' + 'bop_data'
    path += '/' + config.class_type.value
    path += '/' + file_name
    return path


def get_param_path(config, file_name):
    path = config.fs.value + '/' + config.db.value + '/' + 'param'

    path += '/' + config.dt.value
    if config.dt is DataType.gene:

        path += '/' + config.scenario.value
        if config.scenario is Scenario.validation:

            path += '/' + config.validation.value + \
                    '/' + config.validation_method.value
            if config.validation is Validation.simple:
                path += ''

        elif config.scenario is Scenario.approach:
            path += '/' + config.approach.value + \
                    '/' + config.approach_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value

            if config.approach is Approach.top:
                path += '/' + config.approach_gd.value
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                elif config.approach_gd is GeneDataType.from_bop:
                    path += ''
                else:
                    path += '/' + config.geo_type.value
            elif config.approach is Approach.clustering:
                path += '/' + config.approach_gd.value
                if config.approach_gd is GeneDataType.from_cpg:
                    path += ''
                elif config.approach_gd is GeneDataType.from_bop:
                    path += ''
                else:
                    path += '/' + config.geo_type.value

    elif config.dt is DataType.cpg:
        path += '/' + config.scenario.value
        if config.scenario is Scenario.validation:

            path += '/' + config.validation.value + \
                    '/' + config.validation_method.value + \
                    '/' + config.gender.value + \
                    '/' + config.disease.value


        elif config.scenario is Scenario.approach:

            path += '/' + config.scenario.value + \
                    '/' + config.approach.value + \
                    '/' + config.approach_method.value + \
                    '/' + config.gender.value  + \
                    '/' + config.disease.value

    path += '/' + file_name
    return path