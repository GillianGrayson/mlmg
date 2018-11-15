from enum import Enum


class CellPop(Enum):
    sample_id = 'SampleID'
    plasma_blast = 'PlasmaBlast'
    cd8_p = 'CD8pCD28nCD45RAn'
    cd8_naive = 'CD8naive'
    cd4_naive = 'CD4naive'
    cd8_t = 'CD8T'
    cd4_t = 'CD4T'
    nk = 'NK'
    b_cell = 'Bcell'
    mono = 'Mono'
    gran = 'Gran'