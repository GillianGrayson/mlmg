from enum import Enum


class Annotation(Enum):
    cpg = 'ID_REF'
    chr = 'CHR'
    map_info = 'MAPINFO'
    Probe_SNPs = 'Probe_SNPs'
    Probe_SNPs_10 = 'Probe_SNPs_10'
    gene = 'UCSC_REFGENE_NAME'
    class_type = 'Class'
    geo = 'RELATION_TO_UCSC_CPG_ISLAND'
    bop = 'BOP'
    cross_reactive = 'CROSS_R'