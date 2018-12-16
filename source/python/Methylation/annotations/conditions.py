from config.config import *


def cpg_name_condition(config, annotation):
    cpg = annotation[Annotation.cpg.value]

    match = False
    if len(cpg) > 0:
        match = True

    return match


def dna_region_condition(config, annotation):
    dna_region = config.dna_region
    gene = annotation[Annotation.gene.value]

    match = False
    if dna_region is DNARegionType.any:
        match = True
    elif dna_region is DNARegionType.genic:
        if len(gene) > 0:
            match = True
    elif dna_region is DNARegionType.non_genic:
        if len(gene) == 0:
            match = True

    return match

def snp_condition(config, annotation):
    data_type = config.data_type
    snp = config.snp
    curr_snp = annotation[Annotation.Probe_SNPs.value]
    curr_snp_10 = annotation[Annotation.Probe_SNPs_10.value]

    match = True
    if data_type is DataType.cpg:
        if snp is SNPType.snp_excluded or snp is SNPType.snp_excluded_weak or SNPType.snp_cluster:
            if curr_snp != '' or curr_snp_10 != '':
                match = False
    elif data_type is DataType.gene or data_type is DataType.bop:
        if snp is SNPType.snp_excluded_weak:
            if curr_snp != '' or curr_snp_10 != '':
                match = False

    return match

def cross_reactive_condition(config, annotation):
    data_type = config.data_type
    cross_reactive = config.cross_reactive
    cr = annotation[Annotation.cross_reactive.value]

    match = True
    if data_type is DataType.cpg:
        if cross_reactive is CrossReactiveType.cross_reactive_excluded \
                or cross_reactive is CrossReactiveType.cross_reactive_excluded_weak:
            if cr == 1:
                match = False
    elif data_type is DataType.gene or data_type is DataType.bop:
        if cross_reactive is CrossReactiveType.cross_reactive_excluded_weak:
            if cr == 1:
                match = False

    return match

def chromosome_condition(config, annotation):
    chromosome_type = config.chromosome_type
    chr = annotation[Annotation.chr.value]

    match = False
    if chromosome_type is ChromosomeType.all:
        match = True
    elif chromosome_type is ChromosomeType.non_gender:
        if chr != 'X' and chr != 'Y':
            match = True
    elif chromosome_type is ChromosomeType.x:
        if chr == 'X':
            match = True
    elif chromosome_type is ChromosomeType.y:
        if chr == 'Y':
            match = True

    return match

def class_type_condition(config, annotation):
    class_type = config.class_type
    ct = annotation[Annotation.class_type.value]

    match = False
    if class_type is ClassType.any:
        match = True
    elif class_type is ClassType.class_ab:
        classes_values = [ClassType.class_a.value,
                          ClassType.class_b.value]
        if ct in classes_values:
            match = True
    else:
        if ct == class_type.value:
            match = True

    return match

def geo_type_condition(config, annotation):
    geo_type = config.geo_type
    geo = annotation[Annotation.geo.value]

    target_geo = []
    if geo_type is GeoType.islands:
        target_geo.append('Island')
    elif geo_type is GeoType.shores:
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')
    elif geo_type is GeoType.shores_s:
        target_geo.append('S_Shore')
    elif geo_type is GeoType.shores_n:
        target_geo.append('N_Shore')
    elif geo_type is GeoType.islands_shores:
        target_geo.append('Island')
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')

    match = False
    if geo_type is GeoType.any:
        match = True
    else:
        if geo in target_geo:
            match = True

    return match

def is_included(config, annotation):
    excluded = config.excluded
    cpg = annotation[Annotation.cpg.value]

    match = True
    if cpg in excluded:
        match = False

    return match