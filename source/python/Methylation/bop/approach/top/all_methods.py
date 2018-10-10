from config.config import *
from bop.approach.top.manova.top import save_top_manova
from config.method import Method


def top_proc(config, attributes_types, attribute_target):
    if config.method is Method.manova:
        save_top_manova(config, attributes_types, attribute_target)


data_base = DataBase.GSE40279
data_type = DataType.bop

chromosome_type = ChromosomeTypes.non_gender

class_types = [ClassType.class_ab]

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

scenario = Scenario.approach
approach = Approach.top
methods = [Method.manova]

is_clustering = False

attributes_types = [Attribute.age,
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
attribute_target = Attribute.age

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

                is_clustering=is_clustering
            )

            top_proc(config, attributes_types, attribute_target)


