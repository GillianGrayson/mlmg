from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from cpg.approach.top.enet.top import save_top_enet
from cpg.approach.top.linreg.top import save_top_linreg
from cpg.approach.top.linreg_ols.top import save_top_linreg_ols
from cpg.approach.top.linreg_ols_wo_outliers.top import save_top_linreg_ols_wo_outliers
from cpg.approach.top.anova.top import save_top_anova
from cpg.approach.top.spearman.top import save_top_spearman
from cpg.approach.top.moment.top import save_top_moment


def top_proc(config):
    if config.method is Method.enet:
        save_top_enet(config)
    elif config.method is Method.linreg:
        save_top_linreg(config)
    elif config.method is Method.anova:
        save_top_anova(config)
    elif config.method is Method.spearman:
        save_top_spearman(config)
    elif config.method is Method.moment:
        save_top_moment(config)
    elif config.method is Method.linreg_ols:
        save_top_linreg_ols(config)
    elif config.method is Method.linreg_ols_wo_outliers:
        save_top_linreg_ols_wo_outliers(config)

data_base = DataBase.GSE40279
data_type = DataType.cpg

cross_reactive = CrossReactiveType.cross_reactive_included
snp = SNPType.snp_included
chromosome_type = ChromosomeType.non_gender

dna_region = DNARegionType.genic

disease = Disease.any
genders = [Gender.F, Gender.M]

scenario = Scenario.approach
approach = Approach.top
methods = [
    Method.linreg_ols,
]

is_clustering = False

for method in methods:
    print(method.value)
    for gender in genders:
        print('\t' + gender.value)

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

            is_clustering=is_clustering
        )

        top_proc(config)
