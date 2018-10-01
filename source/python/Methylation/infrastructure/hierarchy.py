from config.config import *
import os
from config.method import Method

def create_fs(categories, init_path):
    curr_path = init_path
    for category in categories.reverse():
        for instance in category:

            result_path = class_path + '/' + 'result'
            if not os.path.exists(result_path):
                os.makedirs(result_path)

            # Approach
            scenario_path = result_path + '/' + Scenario.approach.value
            if not os.path.exists(scenario_path):
                os.makedirs(scenario_path)

            for approach in approaches:
                approach_path = scenario_path + '/' + approach.value
                if not os.path.exists(approach_path):
                    os.makedirs(approach_path)

                for method in approach_methods:
                    method_path = approach_path + '/' + method.value
                    if not os.path.exists(method_path):
                        os.makedirs(method_path)

                    for gender in genders:
                        gender_path = method_path + '/' + gender.value
                        if not os.path.exists(gender_path):
                            os.makedirs(gender_path)

                        for disease in diseases:
                            disease_path = gender_path + '/' + disease.value
                            if not os.path.exists(disease_path):
                                os.makedirs(disease_path)

            # Validation
            scenario_path = result_path + '/' + Scenario.validation.value
            if not os.path.exists(scenario_path):
                os.makedirs(scenario_path)


def create_bop_data_fs(config):
    path = get_path(config, '') + 'bop_data'
    if not os.path.exists(path):
        os.makedirs(path)

    dna_regions = [x.value for x in DNARegion] # Level 0
    chromosome_types = [x.value for x in ChromosomeTypes] # Level 1
    class_types = [x.value for x in ClassType] # Level 2

    dna_region_file = open(path + '/' + 'dna_region.txt', 'w')
    for dna_region in dna_regions:
        dna_region_file.write(dna_region + '\n')
        dna_region_path = path + '/' + dna_region
        if not os.path.exists(dna_region_path):
            os.makedirs(dna_region_path)

        chromosome_file = open(dna_region_path + '/' + 'chromosome_types.txt', 'w')
        for chromosome_type in chromosome_types:
            chromosome_file.write(chromosome_type + '\n')
            chromosome_path = dna_region_path + '/' + chromosome_type
            if not os.path.exists(chromosome_path):
                os.makedirs(chromosome_path)

            class_file = open(chromosome_path + '/' + 'class_types.txt', 'w')
            for class_type in class_types:
                class_file.write(class_type + '\n')
                class_path = chromosome_path + '/' + class_type
                if not os.path.exists(class_path):
                    os.makedirs(class_path)

            class_file.close()

        chromosome_file.close()

    dna_region_file.close()


def create_gene_data_fs(config):
    path = get_path(config, '') + 'gene_data'
    if not os.path.exists(path):
        os.makedirs(path)

    dna_regions = [x.value for x in DNARegion] # Level 0
    chromosome_types = [x.value for x in ChromosomeTypes] # Level 1
    gene_data_types = [x.value for x in GeneDataType] # Level 2
    geo_types = [x.value for x in GeoType] # Level 3

    dna_region_file = open(path + '/' + 'dna_region.txt', 'w')
    for dna_region in dna_regions:
        dna_region_file.write(dna_region + '\n')
        dna_region_path = path + '/' + dna_region
        if not os.path.exists(dna_region_path):
            os.makedirs(dna_region_path)

        chromosome_file = open(dna_region_path + '/' + 'chromosome_types.txt', 'w')
        for chromosome_type in chromosome_types:
            chromosome_file.write(chromosome_type + '\n')
            chromosome_path = dna_region_path + '/' + chromosome_type
            if not os.path.exists(chromosome_path):
                os.makedirs(chromosome_path)

            geo_file = open(chromosome_path + '/' + 'geo_types.txt', 'w')
            for geo_type in geo_types:
                geo_file.write(geo_type + '\n')
                geo_path = chromosome_path + '/' + geo_type
                if not os.path.exists(geo_path):
                    os.makedirs(geo_path)

                gene_data_file = open(geo_path + '/' + 'gene_data_types.txt', 'w')
                for gene_data_type in gene_data_types:
                    gene_data_file.write(gene_data_type + '\n')
                    gene_data_path = geo_path + '/' + gene_data_type
                    if not os.path.exists(gene_data_path):
                        os.makedirs(gene_data_path)

                gene_data_file.close()

            geo_file.close()

        chromosome_file.close()

    dna_region_file.close()


def create_result_bop_fs(config):
    db_path = get_path(config, '') + 'result' + '/' + DataType.bop.value

    # Level 0
    scenarios = [Scenario.approach.value,
                 Scenario.validation.value]

    # Level 1
    approaches = [Approach.top.value]
    validations = [Validation.simple.value]

    # Level 2
    approach_methods = [Method.manova.value]
    validation_methods = [Method.match.value]

    # Level 3
    genders = [Gender.any.value,
               Gender.F.value,
               Gender.M.value]

    # Level 4
    diseases = [Disease.any.value,
                Disease.healthy.value,
                Disease.down_syndrome.value]




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
create_bop_data_fs(config)
create_gene_data_fs(config)
