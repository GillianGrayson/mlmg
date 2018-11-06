from config.config import *
from bop.approach.top.manova.top import save_top_manova
from config.method import Method


def top_proc(config):
    if config.method is Method.manova:
        save_top_manova(config)


data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_type = DataType.bop

chromosome_type = ChromosomeTypes.non_gender

class_types = [ClassType.class_ab]

disease = Disease.any
genders = [Gender.any]

scenario = Scenario.approach
approach = Approach.top
methods = [Method.manova]

is_clustering = False

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

for data_base in data_bases:
    for method in methods:
        print(method.value)
        for gender in genders:
            print('\t' + gender.value)
            for class_type in class_types:
                print('\t\t' + class_type.value)

                config = Config(
                    data_base=data_base,
                    data_type=data_type,

                    chromosome_type=chromosome_type,

                    class_type=class_type,

                    disease=disease,
                    gender=gender,

                    scenario=scenario,
                    approach=approach,
                    method=method,

                    attributes_types=attributes_types,
                    attribute_target=attribute_target,

                    is_clustering=is_clustering
                )

                top_proc(config)
