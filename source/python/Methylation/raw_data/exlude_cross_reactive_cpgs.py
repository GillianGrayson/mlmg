from config.config import *
from infrastructure.load.top import *
from config.method import *
from config.types import *
import pandas as pd
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from annotations.bop import get_dict_bop_genes
from infrastructure.path.path import get_origin_path
import math


def exclude_cpgs(config, data_bases, methods):

    fn = get_origin_path(config, 'cross_reactive_cpgs.txt')
    file = open(fn)
    cpgs_cr = file.read().splitlines()
    file.close()


    for data_base in data_bases:
        config.data_base = data_base

        for method_id in range(0, len(methods)):
            method = methods[method_id]

            suffix = ''
            if method == Method.manova:
                target_str = config.attribute_target[0].value
                for target_id in range(1, len(config.attribute_target)):
                    if isinstance(config.attribute_target[target_id], tuple):
                        target_str += '_' + '_x_'.join([x.value for x in config.attribute_target[target_id]])
                    else:
                        target_str += '_' + config.attribute_target[target_id].value

                types_str = '_'.join([x.value for x in config.attributes_types])
                suffix = '_target(' + target_str + ')_exog(' + types_str + ')'

            fn = 'method(' + method.value + ')' + suffix + '.xlsx'
            keys = get_method_order_metrics(method)
            top_dict = load_top_dict(config, keys, fn=fn)

            save_dict = {}
            for key in keys:
                save_dict[key] = []

            for cpg_id in range(0, len(top_dict[keys[0]])):
                cpg = top_dict[keys[0]][cpg_id]
                if cpg not in cpgs_cr:
                    for key in keys:
                        save_dict[key].append(top_dict[key][cpg_id])

            fn = 'method(' + method.value + ')' + suffix + '_wo_cross_reactive.xlsx'
            fn = get_result_path(config, fn)
            df = pd.DataFrame(save_dict)
            writer = pd.ExcelWriter(fn, engine='xlsxwriter')
            df.to_excel(writer, index=False)
            writer.save()


target_data_bases = [DataBase.GSE40279, DataBase.GSE87571]
target_methods = [Method.linreg_ols]

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

attributes_types = [Attribute.gender, Attribute.age]
attribute_target = [Attribute.gender, Attribute.age, (Attribute.gender, Attribute.age)]

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

exclude_cpgs(config, target_data_bases, target_methods)
