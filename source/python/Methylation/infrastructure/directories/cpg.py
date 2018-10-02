from config.types import *
from config.method import *


def get_dna_regions(config):
    dna_regions = [DNARegion.any.value, DNARegion.genic.value]
    return dna_regions