from config.config import *
from infrastructure.load.top import *
from cpg.approach.top.enet.top import save_top_enet
from cpg.approach.top.linreg.top import save_top_linreg
from cpg.approach.top.anova.top import save_top_anova
from cpg.approach.top.spearman.top import save_top_spearman


def top_proc(config):
    if config.method is Method.enet:
        save_top_enet(config)
    elif config.method is Method.linreg:
        save_top_linreg(config)
    elif config.method is Method.anova:
        save_top_anova(config)
    elif config.method is Method.spearman:
        save_top_spearman(config)


data_base = DataBase.GSE40279
data_type = DataType.cpg

chromosome_type = ChromosomeTypes.non_gender

dna_region = DNARegion.genic

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

scenario = Scenario.approach
approach = Approach.top
methods = [
    Method.linreg,
]

is_clustering = False

for method in methods:
    print(method.value)
    for gender in genders:
        print('\t' + gender.value)

        config = Config(
            data_base=data_base,
            data_type=data_type,

            chromosome_type=chromosome_type,

            dna_region=dna_region,

            disease=disease,
            gender=gender,

            scenario=scenario,
            approach=approach,
            method=method,

            is_clustering=is_clustering
        )

        top_proc(config)
