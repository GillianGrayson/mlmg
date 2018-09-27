from config.config import *
import os

from config.method import Method


def create_hierarchy(config):

    db_path = get_path(config, '') + 'result'
    print(db_path)

    data_types = [DataType.bop.value, DataType.gene.value, DataType.cpg.value]

    scenarios = [Scenario.approach.value, Scenario.validation.value]

    validations = [Validation.simple.value]
    approaches = [Approach.top.value, Approach.clustering.value, Approach.statistics.value]

    val_methods = [Method.linreg.value, Method.linreg_mult.value, Method.match.value]
    apr_top_methods = [Method.anova.value,
                       Method.enet.value,
                       Method.linreg.value,
                       Method.linreg_with_rejection.value,
                       Method.linreg_bend.value,
                       Method.linreg_dispersion.value,
                       Method.linreg_variance.value,
                       Method.manova.value,
                       Method.spearman.value]
    apr_bend_methods = [ Method.linreg.value]
    apr_clustering_methods = []

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

                        for method in apr_clustering_methods:

                            method_path = approach_path + '/' + method
                            if not os.path.exists(method_path):
                                os.makedirs(method_path)

                    elif approach == Approach.statistics.value:

                        for gender in genders:

                            gender_path = approach_path + '/' + gender
                            if not os.path.exists(gender_path):
                                os.makedirs(gender_path)

                            for disease in diseases:

                                disease_path = gender_path + '/' + disease
                                if not os.path.exists(disease_path):
                                    os.makedirs(disease_path)

                                for geo in geos:

                                    geo_path = disease_path + '/' + geo
                                    if not os.path.exists(geo_path):
                                        os.makedirs(geo_path)

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

                            for approach in approaches:

                                approach_path = method_path + '/' + approach
                                if not os.path.exists(approach_path):
                                    os.makedirs(approach_path)

                                if approach == Approach.top.value:

                                    for method in apr_top_methods:

                                        method_apr_path = approach_path + '/' + method
                                        if not os.path.exists(method_apr_path):
                                            os.makedirs(method_apr_path)

                                        for gender in genders:

                                            gender_path = method_apr_path + '/' + gender
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

                                                    for geo in geos:

                                                        geo_path = disease_path + '/' + geo
                                                        if not os.path.exists(geo_path):
                                                            os.makedirs(geo_path)

                                                        for gd_o in gds:

                                                            gd_path_order = geo_path + '/' + gd_o
                                                            if not os.path.exists(gd_path_order):
                                                                os.makedirs(gd_path_order)

                                                            for gd_v in gds:

                                                                gd_path_vals = gd_path_order + '/' + gd_v
                                                                if not os.path.exists(gd_path_vals):
                                                                    os.makedirs(gd_path_vals)


config = Config(db=DataBaseType.GSE87571,
                read_only=True)


create_hierarchy(config)
