from enum import Enum

class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    anova = 'anova'
    spearman = 'spearman'
    manova = 'manova'
    clustering = 'clustering'

class Validation(Enum):
    linreg_mult = 'validation_linreg_mult'
    linreg = 'validation_linreg'

class ClusteringType(Enum):
    k_means = 'k_means'
