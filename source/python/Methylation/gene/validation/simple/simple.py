from config.config import *
from gene.validation.simple.linreg_mult.simple import save_simple_linreg_mult
from gene.validation.simple.linreg.simple import save_simple_linreg


def simple_proc(config, num_top=100, num_bootstraps=100):
    if config.validation_method is Method.linreg_mult:
        save_simple_linreg_mult(config, num_top=num_top, num_bootstrap_runs=num_bootstraps)
    elif config.validation_method is Method.linreg:
        save_simple_linreg(config, num_top=num_top)


num_bootstraps = 100
genders = [Gender.F, Gender.any, Gender.M]
tops = range(5, 1001, 5)

for gender in genders:

    config = Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.validation,
        approach=Approach.top,
        validation=Validation.simple,
        approach_method=Method.enet,
        validation_method=Method.linreg_mult,
        approach_gd=GeneDataType.mean,
        validation_gd=GeneDataType.mean,
        geo=GeoType.islands_shores,
        gender=gender
    )

    for num_top in tops:

        simple_proc(config, num_top=num_top, num_bootstraps=num_bootstraps)
