from config.types import *

def regular_condition(config, annotation):
    geo_type = config.geo_type
    dna_region = config.dna_region

    cpg = annotation['ID_REF']
    gene = annotation['UCSC_REFGENE_NAME']
    chr = annotation['CHR']
    geo = annotation['RELATION_TO_UCSC_CPG_ISLAND']

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

    cpg = annotations['ID_REF']
    gene = annotations['UCSC_REFGENE_NAME']
    chr = annotations['CHR']
    geo = annotations['RELATION_TO_UCSC_CPG_ISLAND']
    map_info = annotations['MAPINFO']

    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_geo = geo[i].rstrip()
        curr_map_info = map_info[i].rstrip()

        annotation = {}
        annotation['ID_REF'] = curr_cpg
        annotation['UCSC_REFGENE_NAME'] = curr_gene
        annotation['CHR'] = curr_chr
        annotation['RELATION_TO_UCSC_CPG_ISLAND'] = curr_geo
        annotation['MAPINFO'] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            all_genes = list(set(curr_gene.split(';')))
            dict_cpg_gene[curr_cpg] = all_genes

    return dict_cpg_gene

def get_dict_gene_cpg(config):
    annotations = config.annotations

    cpg = annotations['ID_REF']
    gene = annotations['UCSC_REFGENE_NAME']
    chr = annotations['CHR']
    geo = annotations['RELATION_TO_UCSC_CPG_ISLAND']
    map_info = annotations['MAPINFO']

    dict_gene_cpg = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_geo = geo[i].rstrip()
        curr_map_info = map_info[i].rstrip()

        annotation = {}
        annotation['ID_REF'] = curr_cpg
        annotation['UCSC_REFGENE_NAME'] = curr_gene
        annotation['CHR'] = curr_chr
        annotation['RELATION_TO_UCSC_CPG_ISLAND'] = curr_geo
        annotation['MAPINFO'] = curr_map_info

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

    cpg = annotations['ID_REF']
    gene = annotations['UCSC_REFGENE_NAME']
    chr = annotations['CHR']
    geo = annotations['RELATION_TO_UCSC_CPG_ISLAND']
    map_info = annotations['MAPINFO']

    dict_cpg_map = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_geo = geo[i].rstrip()
        curr_map_info = map_info[i].rstrip()

        annotation = {}
        annotation['ID_REF'] = curr_cpg
        annotation['UCSC_REFGENE_NAME'] = curr_gene
        annotation['CHR'] = curr_chr
        annotation['RELATION_TO_UCSC_CPG_ISLAND'] = curr_geo
        annotation['MAPINFO'] = curr_map_info

        is_passed = regular_condition(config, annotation)

        if is_passed:
            dict_cpg_map[curr_cpg] = curr_map_info

    return dict_cpg_map
