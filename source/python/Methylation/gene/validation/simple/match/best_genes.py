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
genders = [Gender.any, Gender.M, Gender.F]
approach_gd = GeneDataType.mean
geos = [GeoType.islands, GeoType.islands_shores, GeoType.any]

for gt in genders:
    print('gender: ' + gt.value)
    for geo in geos:
        print('\t' + 'geo: ' + geo.value)

        config = Config(
            db=db,
            dt=dt,
            approach=approach,
            scenario=scenario,
            approach_method=Method.linreg,
            gt=gt,
            approach_gd=approach_gd,
            geo=geo
        )

        fn = 'claudio2015.txt'
        gene_top = load_top_gene_names_by_article(config, fn)
        gene_top_dict = {}
        gene_top_dict['claudio2015'] = gene_top

        for method in approach_methods:
            config.approach_method = method
            gene_top = load_top_gene_names(config, num_top)[0:num_top]
            gene_top_dict[method.value] = gene_top

        genes_common_dict = {}
        for s in gene_top_dict:
            genes = gene_top_dict.get(s)
            for gene in genes:
                if gene in genes_common_dict:
                    genes_common_dict[gene] += 1
                else:
                    genes_common_dict[gene] = 1

        genes = list(genes_common_dict.keys())
        counts = list(genes_common_dict.values())
        order = np.argsort(list(map(abs, counts)))[::-1]
        genes_sorted = list(np.array(genes)[order])
        counts_sorted = list(np.array(counts)[order])

        fn = get_result_path(config, 'top.txt')
        save_features(fn, [genes_sorted, counts_sorted])

        config.scenario = Scenario.validation
        config.validation = Validation.simple
        config.validation_method = Method.match
        config.validation_gd = config.approach_gd
        fn = 'best_genes.txt'
        fn = get_result_path(config, fn)
        save_features(fn, [genes_sorted, counts_sorted])