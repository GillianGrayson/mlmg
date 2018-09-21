from config.types import *
import numpy as np
from annotations.regular import get_dict_cpg_gene


def bop_condition(config, annotation):
    cpg_classes = config.cpg_classes
    dna_region = config.dna_region

    cpg = annotation[Annotation.cpg.value]
    chr = annotation[Annotation.chr.value]
    gene = annotation[Annotation.gene.value]
    ct = annotation[Annotation.class_type.value]

    classes_values = [x.values for x in cpg_classes]

    is_class_match = False
    if ct in classes_values:
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

    cpg = annotations[Annotation.cpg.value]
    gene = annotations[Annotation.gene.value]
    chr = annotations[Annotation.chr.value]
    geo = annotations[Annotation.geo.value]
    map_info = annotations[Annotation.map_info.value]
    bop = annotations[Annotation.bop.value]
    class_type = annotations[Annotation.class_type.value]

    dict_bop_cpgs = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i]
        curr_gene = gene[i]
        curr_chr = chr[i]
        curr_geo = geo[i]
        curr_map_info = map_info[i]
        curr_bop = bop[i]
        curr_class_type = class_type[i]

        annotation = {}
        annotation[Annotation.cpg.value] = curr_cpg
        annotation[Annotation.gene.value] = curr_gene
        annotation[Annotation.chr.value] = curr_chr
        annotation[Annotation.geo.value] = curr_geo
        annotation[Annotation.map_info.value] = curr_map_info
        annotation[Annotation.bop.value] = curr_bop
        annotation[Annotation.class_type.value] = curr_class_type

        is_passed = bop_condition(config, annotation)

        if is_passed:
            if len(curr_bop) > 0:
                if curr_bop in dict_bop_cpgs:
                    dict_bop_cpgs[curr_bop].append(curr_cpg)
                else:
                    dict_bop_cpgs[curr_bop] = [curr_cpg]

    # Sorting cpgs in bops by map_info
    num_bops = 0
    for bop in dict_bop_cpgs:
        cpgs = dict_bop_cpgs.get(bop)
        map_infos = []
        for curr_cpg in cpgs:
            cpg_index = cpg.index(curr_cpg)
            curr_map_info = map_info[cpg_index]
            map_infos.append(curr_map_info)
        order = np.argsort(map_infos)
        cpgs_sorted = list(np.array(cpgs)[order])
        dict_bop_cpgs[bop] = cpgs_sorted
        num_bops += 1
        if num_bops % config.print_rate == 0:
            print('num_bops: ' + str(num_bops))

    return dict_bop_cpgs

def get_dict_bop_genes(config, dict_bop_cpgs):
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_bop_genes = {}
    for bop in dict_bop_cpgs:
        cpgs = dict_bop_cpgs.get(bop)
        genes = []
        for curr_cpg in cpgs:
            curr_genes = dict_cpg_gene.get(curr_cpg)
            genes += curr_genes
        dict_bop_genes[bop] = list(set(genes))
    return dict_bop_genes