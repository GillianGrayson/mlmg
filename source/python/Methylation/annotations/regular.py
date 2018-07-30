from config.types import *

def regular_condition(config, annotation):
    geo_type = config.geo_type
    dna_region = config.dna_region

    cpg = annotation[Annotation.cpg.value]
    gene = annotation[Annotation.gene.value]
    chr = annotation[Annotation.chr.value]
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

    is_target_geo = False
    if len(target_geo) == 0:
        is_target_geo = True
    else:
        if geo in target_geo:
            is_target_geo = True

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
    if len(cpg) > 2:
        if cpg[0:2] == 'cg' or cpg[0:2] == 'rs' or cpg[0:2] == 'ch':
            if chr != 'X' and chr != 'Y':
                if is_target_geo:
                    if is_dna_region_match:
                        result = True

    return result

def get_dict_cpg_gene(config):
    annotations = config.annotations

    cpg = annotations[Annotation.cpg.value]
    gene = annotations[Annotation.gene.value]
    chr = annotations[Annotation.chr.value]
    geo = annotations[Annotation.geo.value]
    map_info = annotations[Annotation.map_info.value]

    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i]
        curr_gene = gene[i]
        curr_chr = chr[i]
        curr_geo = geo[i]
        curr_map_info = map_info[i]

        annotation = {}
        annotation[Annotation.cpg.value] = curr_cpg
        annotation[Annotation.gene.value] = curr_gene
        annotation[Annotation.chr.value] = curr_chr
        annotation[Annotation.geo.value] = curr_geo
        annotation[Annotation.map_info.value] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            all_genes = list(set(curr_gene.split(';')))
            dict_cpg_gene[curr_cpg] = all_genes

    return dict_cpg_gene

def get_dict_gene_chr(config):
    annotations = config.annotations

    cpg = annotations[Annotation.cpg.value]
    gene = annotations[Annotation.gene.value]
    chr = annotations[Annotation.chr.value]
    geo = annotations[Annotation.geo.value]
    map_info = annotations[Annotation.map_info.value]

    dict_gene_chr = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i]
        curr_gene = gene[i]
        curr_chr = chr[i]
        curr_geo = geo[i]
        curr_map_info = map_info[i]

        annotation = {}
        annotation[Annotation.cpg.value] = curr_cpg
        annotation[Annotation.gene.value] = curr_gene
        annotation[Annotation.chr.value] = curr_chr
        annotation[Annotation.geo.value] = curr_geo
        annotation[Annotation.map_info.value] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            all_genes = list(set(curr_gene.split(';')))
            for g in all_genes:
                if g in dict_gene_chr:
                    if curr_chr not in dict_gene_chr[g]:
                        dict_gene_chr[g].append(curr_chr)
                else:
                    dict_gene_chr[g] = [curr_chr]

    return dict_gene_chr


def get_dict_gene_cpg(config):
    annotations = config.annotations

    cpg = annotations[Annotation.cpg.value]
    gene = annotations[Annotation.gene.value]
    chr = annotations[Annotation.chr.value]
    geo = annotations[Annotation.geo.value]
    map_info = annotations[Annotation.map_info.value]

    dict_gene_cpg = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i]
        curr_gene = gene[i]
        curr_chr = chr[i]
        curr_geo = geo[i]
        curr_map_info = map_info[i]

        annotation = {}
        annotation[Annotation.cpg.value] = curr_cpg
        annotation[Annotation.gene.value] = curr_gene
        annotation[Annotation.chr.value] = curr_chr
        annotation[Annotation.geo.value] = curr_geo
        annotation[Annotation.map_info.value] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            all_genes = list(set(curr_gene.split(';')))
            for g in all_genes:
                if g in dict_gene_cpg:
                    dict_gene_cpg[g].append(curr_cpg)
                else:
                    dict_gene_cpg[g] = [curr_cpg]

    return dict_gene_cpg

def get_dict_cpg_map_info(config):
    annotations = config.annotations

    cpg = annotations[Annotation.cpg.value]
    gene = annotations[Annotation.gene.value]
    chr = annotations[Annotation.chr.value]
    geo = annotations[Annotation.geo.value]
    map_info = annotations[Annotation.map_info.value]

    dict_cpg_map = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_geo = geo[i].rstrip()
        curr_map_info = map_info[i].rstrip()

        annotation = {}
        annotation[Annotation.cpg.value] = curr_cpg
        annotation[Annotation.gene.value] = curr_gene
        annotation[Annotation.chr.value] = curr_chr
        annotation[Annotation.geo.value] = curr_geo
        annotation[Annotation.map_info.value] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            dict_cpg_map[curr_cpg] = curr_map_info

    return dict_cpg_map
