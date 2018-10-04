from config.config import *
from infrastructure.load.top import *
from cpg.approach.top.enet.top import save_top_enet
from cpg.approach.top.linreg.top import save_top_linreg
from cpg.approach.top.anova.top import save_top_anova
from cpg.approach.top.spearman.top import save_top_spearman


def top_proc(config):
    if config.approach_method is Method.enet:
        save_top_enet(config)
    elif config.approach_method is Method.linreg:
        save_top_linreg(config)
    elif config.approach_method is Method.anova:
        save_top_anova(config)
    elif config.approach_method is Method.spearman:
        save_top_spearman(config)


db = DataBase.GSE40279,
dt = DataType.cpg,
approach = Approach.top,
scenario = Scenario.approach,
approach_methods = [Method.anova, Method.linreg, Method.spearman]
gts = [Gender.any, Gender.M, Gender.F]


for method in approach_methods:
    print('method: ' + method.value)
    for gender in gts:
        print('\t' + 'gender: ' + gender.value)
        config = Config(
            db=DataBase.GSE40279,
            dt=DataType.cpg,
            approach=Approach.top,
            scenario=Scenario.approach,
            approach_method=method,
            gender=gender
        )

        top_proc(config)