from config.config import *
from config.method import *


def get_diseases(config):
    if config.data_base is DataBase.GSE52588:
        diseases = [x.value for x in Disease]
    else:
        diseases = [Disease.any.value]
    return diseases


def get_genders(config):
    genders = [x.value for x in Gender]
    return genders


def get_scenarios(config):
    scenarios = [x.value for x in Scenario]
    return scenarios


def get_approaches(config):
    if config.data_type is DataType.bop:
        approaches = [Approach.top.value]
    elif config.data_type is DataType.cpg:
        approaches = [Approach.top.value]
    elif config.data_type is DataType.gene:
        approaches = [Approach.top.value,
                      Approach.statistics.value]
    else:
        approaches = [x.value for x in Approach]
    return approaches


def get_methods(config):
    if config.data_type is DataType.bop:
        methods = [Method.manova.value,
                   Method.match.value]
    elif config.data_type is DataType.cpg:
        methods = [Method.match.value,
                   Method.anova.value,
                   Method.enet.value,
                   Method.linreg.value,
                   Method.linreg_mult.value,
                   Method.linreg_with_rejection.value,
                   Method.linreg_bend.value,
                   Method.linreg_dispersion.value,
                   Method.linreg_variance.value,
                   Method.spearman.value]
    elif config.data_type is DataType.gene:
        methods = [Method.match.value,
                   Method.anova.value,
                   Method.enet.value,
                   Method.linreg.value,
                   Method.linreg_mult.value,
                   Method.linreg_with_rejection.value,
                   Method.linreg_bend.value,
                   Method.linreg_dispersion.value,
                   Method.linreg_variance.value,
                   Method.spearman.value]
    else:
        methods = [x.value for x in Method]

    return methods


