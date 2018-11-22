from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from cpg.validation.clock.linreg_mult.clock import *


def clock_proc(config_lvl_2, config_lvl_1):
    if config_lvl_2.method is Method.linreg_mult:
        clock_linreg_mult(config_lvl_2, config_lvl_1)

data_bases = [DataBase.GSE87571]
data_type = DataType.cpg

cross_reactives = [CrossReactiveType.cross_reactive_excluded]
snps = [SNPType.snp_excluded]
chromosome_type = ChromosomeType.non_gender

dna_region = DNARegionType.genic

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

scenario = Scenario.validation
approach = Approach.clock
methods = [
    Method.linreg_mult,
]

scenario_lvl_1 = Scenario.approach
approach_lvl_1 = Approach.top
method_lvl_1 = Method.linreg_ols

is_clustering = False

attributes_types = [Attribute.age]

""""
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
"""

attribute_target = [Attribute.age]

for data_base in data_bases:
    print(data_base.value)
    for cross_reactive in cross_reactives:
        print(cross_reactive.value)
        for snp in snps:
            print(snp.value)
            for method in methods:
                print(method.value)
                for gender in genders:
                    print(gender.value)

                    config_lvl_1 = Config(
                        data_base=data_base,
                        data_type=data_type,

                        cross_reactive=cross_reactive,
                        snp=snp,
                        chromosome_type=chromosome_type,

                        dna_region=dna_region,

                        scenario=scenario_lvl_1,
                        approach=approach_lvl_1,
                        method=method_lvl_1,

                        disease=disease,
                        gender=gender,

                        is_clustering=is_clustering,

                        attributes_types=attributes_types,
                        attribute_target=attribute_target
                    )

                    config_lvl_2 = Config(
                        data_base=data_base,
                        data_type=data_type,

                        cross_reactive=cross_reactive,
                        snp=snp,
                        chromosome_type=chromosome_type,

                        dna_region=dna_region,

                        scenario=scenario,
                        approach=approach,
                        method=method,

                        disease=disease,
                        gender=gender,

                        is_clustering=is_clustering,

                        attributes_types=attributes_types,
                        attribute_target=attribute_target
                    )

                    clock_proc(config_lvl_2, config_lvl_1)
