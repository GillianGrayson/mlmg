from config.config import *
from infrastructure.load.top import *
from config.method import *
from config.types import *
import pandas as pd
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
import math


def genes_intersection(config, data_bases, methods, sort_ids, num_top):
    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)



    gene_lists = []
    for data_base in data_bases:
        config.data_base = data_base

        for method_id in range(0, len(methods)):
            method = methods[method_id]

            fn = 'method(' + method.value + ').xlsx'
            keys = get_method_order_metrics(method)
            top_dict = load_top_dict(config, keys, fn=fn)

            order = np.argsort(top_dict[keys[sort_ids[method_id]]])
            for key in keys:
                top_dict[key] = list(np.array(top_dict[key])[order])

            if config.data_type is DataType.gene:
                gene_lists.append(top_dict[keys[0]][0:num_top])
            elif config.data_type is DataType.cpg:
                config.read_only = False
                dict_cpg_gene = get_dict_cpg_gene(config)

                genes = []
                for name in top_dict[keys[0]]:
                    if len(genes) > num_top:
                        break
                    if name in top_dict[keys[0]]:
                        curr_genes = dict_cpg_gene[name]
                        for gene in curr_genes:
                            if gene not in genes:
                                genes.append(gene)
                gene_lists.append(genes[0:num_top])

    gene_intersection = gene_lists[0]
    for genes in gene_lists:
        gene_intersection = list(set(gene_intersection).intersection(genes))

    gene_intersection.sort()

    for gene in gene_intersection:
        print(gene)
    print(str(len(gene_intersection)))

target_data_bases = [DataBase.GSE87571, DataBase.GSE40279]
target_methods = [Method.linreg_ols]
target_sort_ids = [1]
target_num_top = 750

data_base = DataBase.data_base_versus

data_type = DataType.cpg

chromosome_type = ChromosomeTypes.non_gender

# bop
class_type = ClassType.class_ab
# cpg
dna_region = DNARegion.genic
# gene
geo_type = GeoType.islands_shores
gene_data_type = GeneDataType.mean

scenario = Scenario.validation
approach = Approach.top
method = Method.gender_specific

disease = Disease.any
gender = Gender.versus

is_clustering = False

config = Config(
    read_only=True,

    data_base=data_base,

    data_type=data_type,

    chromosome_type=chromosome_type,

    class_type=class_type,

    dna_region=dna_region,

    geo_type=geo_type,
    gene_data_type=gene_data_type,

    scenario=scenario,
    approach=approach,
    method=method,

    disease=disease,
    gender=gender,

    is_clustering=is_clustering
)

genes_intersection(config, target_data_bases, target_methods, target_sort_ids, target_num_top)
