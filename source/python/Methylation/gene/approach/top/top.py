from config.config import *
from infrastructure.load.top import *
from gene.approach.top.enet.params import save_params_enet
from gene.approach.top.enet.top import save_top_enet
from gene.approach.top.linreg.top import save_top_linreg
from gene.approach.top.anova.top import save_top_anova
from gene.approach.top.spearman.top import save_top_spearman


def top_proc(config):
    if config.approach_method is Method.enet:
        save_top_enet(config)
    elif config.approach_method is Method.linreg:
        save_top_linreg(config)
    elif config.approach_method is Method.anova:
        save_top_anova(config)
    elif config.approach_method is Method.spearman:
        save_top_spearman(config)

config = Config(
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    approach=Approach.top,
    scenario=Scenario.approach,
    approach_method=Method.spearman,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands,
)

top_proc(config)

