from enum import Enum


class CrossReactiveType(Enum):
    cross_reactive_included = 'cross_reactive_included'
    cross_reactive_excluded = 'cross_reactive_excluded'


class SNPType(Enum):
    snp_included = 'snp_included'
    snp_excluded = 'snp_excluded'


class ChromosomeType(Enum):
    all = 'all'
    non_gender = 'non_gender'
    x = 'x'
    y = 'y'