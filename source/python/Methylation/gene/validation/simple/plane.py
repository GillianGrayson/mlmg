from config.config import *
from config.method import Method
from gene.validation.simple.linreg.plane import save_plane_linreg

def plane_proc(config):
    if config.validation_method is Method.linreg:
        save_plane_linreg(config)


config = Config(
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    approach=Approach.top,
    validation=Validation.simple,
    scenario=Scenario.validation,
    approach_method=Method.linreg,
    validation_method=Method.linreg_mult,
    approach_gd=GeneDataType.mean,
    validation_gd=GeneDataType.mean,
    geo=GeoType.any
)

plane_proc(config)