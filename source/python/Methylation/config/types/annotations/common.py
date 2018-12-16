from enum import Enum


class CrossReactiveType(Enum):
    cross_reactive_included = 'cross_reactive_included'
    cross_reactive_excluded = 'cross_reactive_excluded'
    cross_reactive_excluded_weak = 'cross_reactive_excluded_weak'


class SNPType(Enum):
    snp_included = 'snp_included'
    snp_excluded = 'snp_excluded'
    snp_excluded_weak = 'snp_excluded_weak'
    snp_cluster = 'snp_cluster'


class ChromosomeType(Enum):
    all = 'all'
    non_gender = 'non_gender'
    x = 'x'
    y = 'y'