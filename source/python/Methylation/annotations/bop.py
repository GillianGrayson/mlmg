from config.types import *
import numpy as np
from annotations.gene import get_dict_cpg_gene
from annotations.conditions import *
from infrastructure.path.path import *
import os.path
import pickle


def bop_condition(config, annotation):
    match = False
    if cpg_name_condition(config, annotation):
        if dna_region_condition(config, annotation):
            if chromosome_condition(config, annotation):
                if class_type_condition(config, annotation):
                    match = True
    return match


def get_dict_bop_cpgs(config):
    fn_pkl = 'dict_bop_cpgs.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_bop_cpgs = pickle.load(f)
        f.close()
    else:
        dict_bop_cpgs = {}

        annotations = config.annotations

        cpg = annotations[Annotation.cpg.value]
        gene = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        geo = annotations[Annotation.geo.value]
        map_info = annotations[Annotation.map_info.value]
        bop = annotations[Annotation.bop.value]
        class_type = annotations[Annotation.class_type.value]
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

            if bop_condition(config, annotation):
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

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_bop_cpgs, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    return dict_bop_cpgs


def get_dict_bop_genes(config, dict_bop_cpgs):
    fn_pkl = 'dict_bop_genes.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_bop_genes = pickle.load(f)
        f.close()
    else:
        dict_bop_genes = {}
        dict_cpg_gene = get_dict_cpg_gene(config)
        for bop in dict_bop_cpgs:
            cpgs = dict_bop_cpgs.get(bop)
            genes = []
            for curr_cpg in cpgs:
                curr_genes = dict_cpg_gene.get(curr_cpg)
                genes += curr_genes
            dict_bop_genes[bop] = list(set(genes))

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_bop_genes, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    return dict_bop_genes