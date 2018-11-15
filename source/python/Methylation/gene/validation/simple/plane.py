from config.config import *
from config.types.experiments.method import Method
from config.types.annotations import GeneDataType, GeoType
from gene.validation.simple.linreg.plane import save_plane_linreg

def plane_proc(config):
    if config.validation_method is Method.linreg:
        save_plane_linreg(config)


config = Config(
    db=DataBase.GSE40279,
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