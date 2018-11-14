from enum import Enum


class DNARegionType(Enum):
    genic = 'genic'
    non_genic = 'non_genic'
    any = 'any'