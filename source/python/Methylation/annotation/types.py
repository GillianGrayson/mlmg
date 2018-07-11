from enum import Enum

class GeoType(Enum):
    shores = '_shores'
    shores_s = '_shores_s'
    shores_n = '_shores_n'
    islands = '_islands'
    islands_shores = '_islands_shores'
    any = ''

class ClassType(Enum):
    class_a = 'ClassA'
    class_b = 'ClassB'
    class_c = 'ClassC'
    class_d = 'ClassD'
    any = ''

class DNARegion(Enum):
    genic = 'genic'
    non_genic = 'genic'
    any = ''