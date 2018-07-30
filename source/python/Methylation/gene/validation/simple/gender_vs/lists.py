from config.config import *
from infrastructure.load.top import *
from infrastructure.save.features import save_features
import numpy as np

num_top = 100

db = DataBaseType.GSE40279
dt = DataType.gene
approach = Approach.top
scenario = Scenario.approach
approach_methods = [Method.anova, Method.enet, Method.linreg, Method.spearman]
approach_gds = [GeneDataType.mean, GeneDataType.from_cpg]
geos = [GeoType.islands, GeoType.islands_shores, GeoType.any]

for method in approach_methods:
    print('method: ' + method.value)
    for gd in approach_gds:
        print('\t' + 'gd: ' + gd.value)
        for geo in geos:
            print('\t\t' + 'geo: ' + geo.value)
            config = Config(
                db=db,
                dt=dt,
                approach=approach,
                scenario=scenario,
                approach_method=Method.linreg,
                gender=Gender.any,
                approach_gd=gd,
                geo=geo
            )

            config.gender = Gender.M
            genes_m = load_top_gene_names(config, num_top)[0:num_top]

            config.gender = Gender.F
            genes_f = load_top_gene_names(config, num_top)[0:num_top]

            intersection_genes = list(set(genes_m).intersection(genes_f))
            only_m = list(set(genes_m) - set(intersection_genes))
            only_f = list(set(genes_f) - set(intersection_genes))

            gene_top_dict = {}

            config.scenario = Scenario.validation
            config.validation = Validation.simple
            config.validation_method = Method.gender_vs
            config.validation_gd = config.approach_gd
            fn = 'intersection.txt'
            fn = get_result_path(config, fn)
            save_features(fn, [intersection_genes])
            fn = 'only_m.txt'
            fn = get_result_path(config, fn)
            save_features(fn, [only_m])
            fn = 'only_f.txt'
            fn = get_result_path(config, fn)
            save_features(fn, [only_f])