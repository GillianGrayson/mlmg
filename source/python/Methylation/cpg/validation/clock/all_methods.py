from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from cpg.validation.clock.linreg_mult.clock import *
from config.types.auxiliary import *


def clock_proc(config_lvl_2, config_lvl_1):
    if config_lvl_2.method is Method.linreg_mult:
        clock_linreg_mult(config_lvl_2, config_lvl_1)

data_bases = [DataBase.liver]
data_type = DataType.cpg

cross_reactives = [CrossReactiveType.cross_reactive_excluded]
snps = [SNPType.snp_excluded]
chromosome_type = ChromosomeType.non_gender

dna_region = DNARegionType.genic

disease = Disease.any
genders = [Gender.M, Gender.F, Gender.any]

scenario = Scenario.validation
approach = Approach.clock
methods = [
    Method.linreg_mult,
]
method_params = [
    {'exog_type' : ClockExogType.all,
     'exog_num' : 100,
     'exog_num_comb' : 100}
]

scenario_lvl_1 = Scenario.approach
approach_lvl_1 = Approach.top
method_lvl_1 = Method.linreg_ols
lvl_1_method_params = {'outliers_limit': 0.0,
                       'outliers_sigma': 0.0}

is_clustering = False

attributes_types = [Attribute.age]

"""
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
                        attribute_target=attribute_target,

                        method_params=lvl_1_method_params
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
                        attribute_target=attribute_target,

                        method_params = method_params[methods.index(method)]
                    )

                    clock_proc(config_lvl_2, config_lvl_1)
