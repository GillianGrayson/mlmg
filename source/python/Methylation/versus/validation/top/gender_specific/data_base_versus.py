from config.config import *
from infrastructure.load.top import *
from config.method import *
from config.types import *
import pandas as pd
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from annotations.bop import get_dict_bop_genes
import math


def genes_intersection(config, data_bases, methods, sort_ids, sort_directions, num_top):
    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)

    gene_lists = []

    gene_bop = {}
    gene_cpg = {}
    basic_genes = []
    for data_base in data_bases:
        config.data_base = data_base

        for method_id in range(0, len(methods)):
            method = methods[method_id]

            suffix = ''
            if method == Method.manova:
                target_str = config.attribute_target[0].value
                is_cells = False
                for target_id in range(1, len(config.attribute_target)):
                    if isinstance(config.attribute_target[target_id], tuple):
                        tmp = []
                        is_mult_cells = False
                        for x in config.attribute_target[target_id]:
                            if isinstance(x, Attribute):
                                tmp.append(x.value)
                            elif isinstance(x, CellPop):
                                is_mult_cells = True
                        target_str += '_' + '_x_'.join([x for x in tmp])
                        if is_mult_cells:
                            target_str += '_x_cells'
                    elif isinstance(config.attribute_target[target_id], CellPop):
                        if not is_cells:
                            is_cells = True
                            target_str += '_cells'
                    else:
                        target_str += '_' + config.attribute_target[target_id].value

                types_str = config.attributes_types[0].value
                is_cells = False
                for type_id in range(1, len(config.attributes_types)):
                    if isinstance(config.attributes_types[type_id], Attribute):
                        types_str += '_' + config.attributes_types[type_id].value
                    elif isinstance(config.attributes_types[type_id], CellPop):
                        if not is_cells:
                            is_cells = True
                            types_str += '_cells'

                suffix = '_target(' + target_str + ')_exog(' + types_str + ')'

            fn = 'method(' + method.value + ')' + suffix + '.xlsx'
            keys = get_method_order_metrics(method)
            top_dict = load_top_dict(config, keys, fn=fn)

            order = np.argsort(top_dict[keys[sort_ids[method_id]]])
            if sort_directions[method_id] < 0:
                order = order[::-1]

            for key in keys:
                top_dict[key] = list(np.array(top_dict[key])[order])

            if config.data_type is DataType.gene:
                gene_lists.append(top_dict[keys[0]][0:num_top])
                basic_genes = top_dict[keys[0]]
            elif config.data_type is DataType.cpg:
                config.read_only = False
                dict_cpg_gene = get_dict_cpg_gene(config)
                genes = []
                num_passed_names = 0
                for name in top_dict[keys[0]]:
                    num_passed_names += 1
                    if len(genes) > num_top:
                        break
                    if name in top_dict[keys[0]]:
                        curr_genes = dict_cpg_gene[name]
                        for gene in curr_genes:
                            if gene not in genes:
                                gene_cpg[gene] = name
                                genes.append(gene)
                gene_lists.append(genes[0:num_top])
                basic_genes = genes
                print('num_passed_names: ' + str(num_passed_names))
            elif config.data_type is DataType.bop:
                config.read_only = False
                dict_bop_gene = get_dict_bop_genes(config)
                genes = []
                num_passed_names = 0
                for name in top_dict[keys[0]]:
                    num_passed_names += 1
                    if len(genes) > num_top:
                        break
                    if name in top_dict[keys[0]]:
                        curr_genes = dict_bop_gene[name]
                        for gene in curr_genes:
                            if gene not in genes:
                                gene_bop[gene] = name
                                genes.append(gene)
                gene_lists.append(genes[0:num_top])
                basic_genes = genes
                print('num_passed_names: ' + str(num_passed_names))

    gene_intersection = gene_lists[0]
    for genes in gene_lists:
        gene_intersection = list(set(gene_intersection).intersection(genes))

    # Order of intersection genes corresponds to last database
    for gene in basic_genes:
        if gene in gene_intersection:
            if config.data_type is DataType.gene:
                print(gene)
            elif config.data_type is DataType.cpg:
                print(gene + ' ' + gene_cpg[gene])
            elif config.data_type is DataType.bop:
                print(gene + ' ' + gene_bop[gene])

    print(str(len(gene_intersection)))

target_data_bases = [DataBase.GSE40279, DataBase.GSE87571]
target_methods = [Method.linreg_ols]
target_sort_ids = [3, 1]
target_sort_directions = [-1, 1]
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

attributes_types = [Attribute.gender,
                    Attribute.age,
                    CellPop.plasma_blast,
                    CellPop.cd8_p,
                    CellPop.cd4_naive,
                    CellPop.cd8_naive,
                    CellPop.cd8_t,
                    CellPop.cd4_t,
                    CellPop.nk,
                    CellPop.b_cell,
                    CellPop.mono,
                    CellPop.gran]

attribute_target = [Attribute.gender,
                    Attribute.age,
                    (Attribute.gender,
                     Attribute.age)]

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

    attributes_types=attributes_types,
    attribute_target=attribute_target,

    is_clustering=is_clustering
)

genes_intersection(config, target_data_bases, target_methods, target_sort_ids, target_sort_directions, target_num_top)
