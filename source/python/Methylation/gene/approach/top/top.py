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


def top_proc(config):
    if config.approach_method is Method.enet:
        save_params_enet(config)
        save_top_enet(config)
    elif config.approach_method is Method.linreg:
        save_top_linreg(config)
    elif config.approach_method is Method.anova:
        save_top_anova(config)
    elif config.approach_method is Method.spearman:
        save_top_spearman(config)
    elif config.approach_method is Method.linreg_with_rejection:
        save_top_linreg_with_rejection(config)
    elif config.approach_method is Method.linreg_bend:
        save_top_linreg_bend(config)
    elif config.approach_method is Method.linreg_dispersion:
        save_top_linreg_with_dispersion(config)
    elif config.approach_method is Method.linreg_variance:
        save_top_linreg_variance(config)


db = DataBaseType.GSE87571
dt = DataType.gene
approach = Approach.top
scenario = Scenario.approach
approach_methods = [Method.enet, Method.anova, Method.linreg, Method.linreg_variance, Method.linreg_bend, Method.linreg_dispersion, Method.linreg_with_rejection, Method.spearman]
approach_gd = GeneDataType.mean
genders = [Gender.F, Gender.M, Gender.any]
geos = [GeoType.islands_shores]

for method in approach_methods:
    print(method.value)
    for gender in genders:
        print('\t' + gender.value)
        for geo in geos:
            print('\t\t' + geo.value)

            config = Config(
                db=db,
                dt=dt,
                approach=approach,
                scenario=scenario,
                approach_method=method,
                approach_gd=approach_gd,
                geo=geo,
                gender=gender
            )

            top_proc(config)

