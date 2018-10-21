from config.config import *
from infrastructure.load.top import *
from gene.approach.top.enet.params import save_params_enet
from gene.approach.top.enet.top import save_top_enet
from gene.approach.top.linreg.top import save_top_linreg
from gene.approach.top.linreg_with_rejection.top import save_top_linreg_with_rejection
from gene.approach.top.anova.top import save_top_anova
from gene.approach.top.spearman.top import save_top_spearman
from gene.approach.top.linreg_bend.top import save_top_linreg_bend
from gene.approach.top.linreg_dispersion.top import save_top_linreg_with_dispersion
from gene.approach.top.linreg_variance.top import save_top_linreg_variance
from gene.approach.top.moment.top import save_top_moment
from gene.approach.top.linreg_ols.top import save_top_linreg_ols


def top_proc(config):
    if config.method is Method.enet:
        save_params_enet(config)
        save_top_enet(config, )
    elif config.method is Method.linreg:
        save_top_linreg(config)
    elif config.method is Method.anova:
        save_top_anova(config)
    elif config.method is Method.spearman:
        save_top_spearman(config)
    elif config.method is Method.linreg_with_rejection:
        save_top_linreg_with_rejection(config)
    elif config.method is Method.linreg_bend:
        save_top_linreg_bend(config)
    elif config.method is Method.linreg_dispersion:
        save_top_linreg_with_dispersion(config)
    elif config.method is Method.linreg_variance:
        save_top_linreg_variance(config)
    elif config.method is Method.moment:
        save_top_moment(config)
    elif config.method is Method.linreg_ols:
        save_top_linreg_ols(config)

data_base = DataBase.GSE40279
data_type = DataType.gene

chromosome_type = ChromosomeTypes.non_gender

geo_types = [GeoType.islands_shores]
gene_data_type = GeneDataType.mean

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

scenario = Scenario.approach
approach = Approach.top
methods = [
    Method.moment,
]

is_clustering = False

for method in methods:
    print(method.value)
    for gender in genders:
        print('\t' + gender.value)
        for geo_type in geo_types:
            print('\t\t' + geo_type.value)

            config = Config(
                data_base=data_base,
                data_type=data_type,

                chromosome_type=chromosome_type,

                geo_type=geo_type,
                gene_data_type=gene_data_type,

                scenario=scenario,
                approach=approach,
                method=method,

                disease=disease,
                gender=gender,

                is_clustering=is_clustering
            )

            top_proc(config)

