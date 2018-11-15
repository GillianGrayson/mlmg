from enum import Enum


class GeneDataType(Enum):
    mean = 'mean'
    std = 'std'
    mean_der = 'mean_der'
    mean_der_normed = 'mean_der_normed'
    from_cpg = 'from_cpg'
    from_bop = 'from_bop'


class GeoType(Enum):
    shores = 'shores'
    shores_s = 'shores_s'
    shores_n = 'shores_n'
    islands = 'islands'
    islands_shores = 'islands_shores'
    any = 'any'
    from_cpg = 'from_cpg'
    from_bop = 'from_bop'