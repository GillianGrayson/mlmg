from config.config import *
from infrastructure.load.top import *
from config.method import *
from config.types import *
import pandas as pd
from infrastructure.save.features import save_features
import math


def genes_intersection(config, data_bases, method, sort_id, part):
    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)
    gene_lists = []
    for data_base in data_bases:
        config.data_base = data_base

        fn = 'method(' + method.value + ').xlsx'
        keys = get_method_order_metrics(method)
        top_dict = load_top_dict(config, keys, fn=fn)

        part_int = math.floor(len(top_dict[keys[0]]) * part)

        order = np.argsort(top_dict[keys[sort_id]])
        for key in keys:
            top_dict[key] = list(np.array(top_dict[key])[order])

        gene_lists.append(top_dict[keys[0]][0:part_int])


    gene_intersection = gene_lists[0]
    for genes in gene_lists:
        gene_intersection = list(set(gene_intersection).intersection(genes))

    for gene in gene_intersection:
        print(gene)
    print(str(len(gene_intersection)))

target_method = Method.linreg_ols
target_sort_id = 1
target_part = 0.00665335994

data_base = DataBase.data_base_versus
data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_type = DataType.gene

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

genes_intersection(config, data_bases, target_method, target_sort_id, target_part)
