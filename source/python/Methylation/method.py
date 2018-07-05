from enum import Enum

class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    anova = 'anova'
    spearman = 'spearman'

class Validation(Enum):
    linreg_mult = 'val_linreg_mult'
