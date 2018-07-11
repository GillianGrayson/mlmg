from enum import Enum

class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    anova = 'anova'
    spearman = 'spearman'
    manova = 'manova'

class Validation(Enum):
    linreg_mult = 'validation_linreg_mult'
    linreg = 'validation_linreg'
