from annotation.types import *

def bop_condition(config, annotation):
    target_ct = config.class_type
    dna_region = config.dna_region

    cpg = annotation['ID_REF']
    chr = annotation['CHR']
    gene = annotation['UCSC_REFGENE_NAME']
    ct = annotation['Class']

    is_class_match = False
    if ct is ClassType.any or ct is target_ct:
        is_class_match = True

    is_dna_region_match = False
    if dna_region is DNARegion.any:
        is_dna_region_match = True
    elif dna_region is DNARegion.genic:
        if len(gene) > 0:
            is_dna_region_match = True
    elif dna_region is DNARegion.non_genic:
        if len(gene) == 0:
            is_dna_region_match = True

    result = False
    if is_class_match:
        if is_dna_region_match:
            if len(cpg) > 2:
                if cpg[0:2] == 'cg' or cpg[0:2] == 'rs' or cpg[0:2] == 'ch':
                    if chr != 'X' and chr != 'Y':
                        result = True

    return result


def get_dict_bop_cpgs(config):
    annotations = config.annotations

    cpg = annotations['ID_REF']
    gene = annotations['UCSC_REFGENE_NAME']
    chr = annotations['CHR']
    geo = annotations['RELATION_TO_UCSC_CPG_ISLAND']
    map_info = annotations['MAPINFO']
    bop = annotations['BOP']
    class_type = annotations['Class']

    dict_bop_cpgs = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_geo = geo[i].rstrip()
        curr_map_info = map_info[i].rstrip()
        curr_bop = bop[i].rstrip()
        curr_class_type = class_type[i].rstrip()

        annotation = {}
        annotation['ID_REF'] = curr_cpg
        annotation['UCSC_REFGENE_NAME'] = curr_gene
        annotation['CHR'] = curr_chr
        annotation['RELATION_TO_UCSC_CPG_ISLAND'] = curr_geo
        annotation['MAPINFO'] = curr_map_info
        annotation['BOP'] = curr_bop
        annotation['Class'] = curr_class_type

        is_passed = bop_condition(config, annotation)

        if len(curr_bop) > 0:
            if curr_bop in dict_bop_cpgs:
                if is_passed:
                    dict_bop_cpgs[curr_bop].append(curr_cpg)
            else:
                dict_bop_cpgs[curr_bop] = [curr_cpg]

    return dict_bop_cpgs
