from config.config import *
from gene.validation.simple.linreg_mult.simple import save_simple_linreg_mult
from gene.validation.simple.linreg.simple import save_simple_linreg


def simple_proc(config):
    if config.validation_method is Method.linreg_mult:
        save_simple_linreg_mult(config)
    elif config.validation_method is Method.linreg:
        save_simple_linreg(config)


config = Config(
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    approach=Approach.top,
    validation=Validation.simple,
    scenario=Scenario.validation,
    approach_method=Method.enet,
    validation_method=Method.linreg_mult,
    approach_gd=GeneDataType.mean,
    validation_gd=GeneDataType.mean,
    geo=GeoType.islands
)

simple_proc(config)
