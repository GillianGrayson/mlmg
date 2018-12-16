from infrastructure.path.path import *
from annotations.conditions import *
import os.path
import pickle
import numpy as np


def cpg_condition(config, annotation):
    match = False
    if is_included(config, annotation):
        if cross_reactive_condition(config, annotation):
            if snp_condition(config, annotation):
                if chromosome_condition(config, annotation):
                    if cpg_name_condition(config, annotation):
                        if dna_region_condition(config, annotation):
                            match = True
    return match


def get_approved_cpgs(config):
    fn_pkl = 'approved_cpgs.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        approved_cpgs = pickle.load(f)
        f.close()
    else:
        approved_cpgs = []

        annotations = config.annotations

        cpg = annotations[Annotation.cpg.value]
        gene = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        snp = annotations[Annotation.Probe_SNPs.value]
        snp1_10 = annotations[Annotation.Probe_SNPs_10.value]
        cross_reactive = annotations[Annotation.cross_reactive.value]
        for i in range(0, len(cpg)):

            curr_cpg = cpg[i]
            curr_gene = gene[i]
            curr_chr = chr[i]
            curr_snp = snp[i]
            curr_snp_10 = snp1_10[i]
            curr_cross_reactive = cross_reactive[i]

            annotation = {}
            annotation[Annotation.cpg.value] = curr_cpg
            annotation[Annotation.gene.value] = curr_gene
            annotation[Annotation.chr.value] = curr_chr
            annotation[Annotation.Probe_SNPs.value] = curr_snp
            annotation[Annotation.Probe_SNPs_10.value] = curr_snp_10
            annotation[Annotation.cross_reactive.value] = curr_cross_reactive

            if cpg_condition(config, annotation):
                approved_cpgs.append(curr_cpg)

        f = open(fn_pkl, 'wb')
        pickle.dump(approved_cpgs, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    return approved_cpgs
