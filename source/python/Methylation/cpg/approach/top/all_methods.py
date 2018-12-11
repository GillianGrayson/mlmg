from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from cpg.approach.top.enet.top import save_top_enet
from cpg.approach.top.linreg.top import save_top_linreg
from cpg.approach.top.linreg_ols.top import save_top_linreg_ols
from cpg.approach.top.linreg_variance_ols.top import save_top_linreg_variance_ols
from cpg.approach.top.anova.top import save_top_anova
from cpg.approach.top.anova_statsmodels.top import save_top_anova_statsmodels
from cpg.approach.top.spearman.top import save_top_spearman
from cpg.approach.top.moment.top import save_top_moment


def top_proc(config):
    if config.method is Method.enet:
        save_top_enet(config)
    elif config.method is Method.linreg:
        save_top_linreg(config)
    elif config.method is Method.anova:
        save_top_anova(config)
    elif config.method is Method.anova_statsmodels:
        save_top_anova_statsmodels(config)
    elif config.method is Method.spearman:
        save_top_spearman(config)
    elif config.method is Method.moment:
        save_top_moment(config)
    elif config.method is Method.linreg_ols:
        save_top_linreg_ols(config)
    elif config.method is Method.linreg_variance_ols:
        save_top_linreg_variance_ols(config)

data_bases = [DataBase.liver]
data_type = DataType.cpg

cross_reactives = [CrossReactiveType.cross_reactive_excluded]
snps = [SNPType.snp_excluded]
chromosome_type = ChromosomeType.non_gender

dna_region = DNARegionType.genic

disease = Disease.any
genders = [Gender.any, Gender.M, Gender.F]

scenario = Scenario.approach
approach = Approach.top
methods = [
    Method.linreg_ols,
]
method_params = [
    {'outliers_limit': 0.0,
     'outliers_sigma': 0.0}
]

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

                    config = Config(
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

                        method_params=method_params[methods.index(method)]
                    )

                    top_proc(config)
