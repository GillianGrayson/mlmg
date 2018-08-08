from config.config import *
import os


def create_hierarchy(config):

    db_path = get_path(config, '') + 'result'

    data_types = [DataType.bop.value, DataType.gene.value, DataType.cpg.value]

    scenarios = [Scenario.approach.value, Scenario.validation.value]

    validations = [Validation.simple.value]
    approaches = [Approach.top.value, Approach.clustering.value]

    val_methods = [Method.linreg.value, Method.linreg_mult.value, Method.match.value]
    apr_top_methods = [Method.anova.value, Method.enet.value, Method.linreg.value, Method.manova.value, Method.spearman.value]
    apr_cluster_methods = []

    genders = [Gender.any.value, Gender.F.value, Gender.M.value]

    diseases = [Disease.any.value, Disease.healthy.value, Disease.down_syndrome.value]

    gds = [GeneDataType.mean.value, GeneDataType.mean_der.value, GeneDataType.mean_der_normed.value, GeneDataType.std.value, GeneDataType.from_cpg.value, GeneDataType.from_bop.value]

    geos = [GeoType.any.value, GeoType.islands.value, GeoType.islands_shores.value]

    for data_type in data_types:

        data_type_path = db_path + '/' + data_type
        if not os.path.exists(data_type_path):
            os.makedirs(data_type_path)

        for scenario in scenarios:

            scenario_path = data_type_path + '/' + scenario
            if not os.path.exists(scenario_path):
                os.makedirs(scenario_path)

            if scenario == Scenario.approach.value:

                for approach in approaches:

                    approach_path = scenario_path + '/' + approach
                    if not os.path.exists(approach_path):
                        os.makedirs(approach_path)

                    if approach == Approach.top.value:

                        for method in apr_top_methods:

                            method_path = approach_path + '/' + method
                            if not os.path.exists(method_path):
                                os.makedirs(method_path)

                            for gender in genders:

                                gender_path = method_path + '/' + gender
                                if not os.path.exists(gender_path):
                                    os.makedirs(gender_path)

                                for disease in diseases:

                                    disease_path = gender_path + '/' + disease
                                    if not os.path.exists(disease_path):
                                        os.makedirs(disease_path)

                                    if data_type == DataType.bop.value:

                                        pass

                                    elif data_type == DataType.cpg.value:

                                        pass

                                    elif data_type == DataType.gene.value:

                                        for gd in gds:

                                            gd_path = disease_path + '/' + gd
                                            if not os.path.exists(gd_path):
                                                os.makedirs(gd_path)

                                            if gd == GeneDataType.from_cpg.value or gd == GeneDataType.from_bop.value:

                                                pass

                                            else:

                                                for geo in geos:

                                                    geo_path = gd_path + '/' + geo
                                                    if not os.path.exists(geo_path):
                                                        os.makedirs(geo_path)

                    elif approach == Approach.clustering.value:

                        for method in apr_cluster_methods:

                            method_path = approach_path + '/' + method
                            if not os.path.exists(method_path):
                                os.makedirs(method_path)

            elif scenario == Scenario.validation.value:

                for validation in validations:

                    validation_path = scenario_path + '/' + validation
                    if not os.path.exists(validation_path):
                        os.makedirs(validation_path)

                    for method in val_methods:

                        method_path = validation_path + '/' + method
                        if not os.path.exists(method_path):
                            os.makedirs(method_path)

                        if method == Method.match.value:

                            for gender in genders:

                                gender_path = method_path + '/' + gender
                                if not os.path.exists(gender_path):
                                    os.makedirs(gender_path)

                                for disease in diseases:

                                    disease_path = gender_path + '/' + disease
                                    if not os.path.exists(disease_path):
                                        os.makedirs(disease_path)

                                    if data_type == DataType.bop.value:

                                        pass

                                    elif data_type == DataType.cpg.value:

                                        pass

                                    elif data_type == DataType.gene.value:

                                        for gd in gds:

                                            gd_path = disease_path + '/' + gd
                                            if not os.path.exists(gd_path):
                                                os.makedirs(gd_path)

                                            if gd == GeneDataType.from_cpg.value or gd == GeneDataType.from_bop.value:

                                                pass

                                            else:

                                                for geo in geos:

                                                    geo_path = gd_path + '/' + geo
                                                    if not os.path.exists(geo_path):
                                                        os.makedirs(geo_path)

                        elif method == Method.linreg.value:

                            pass

                        elif method == Method.linreg_mult.value:

                            pass

config = Config(db=DataBaseType.GSE61256,
                read_only=True)

create_hierarchy(config)
