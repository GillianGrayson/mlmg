from infrastructure.path.path import *
from annotations.conditions import *
import os.path
import pickle
import numpy as np


def gene_condition(config, annotation):
    match = False
    if chromosome_condition(config, annotation):
        if cpg_name_condition(config, annotation):
            if dna_region_condition(config, annotation):
                if geo_type_condition(config, annotation):
                    match = True
    return match


def get_dict_cpg_gene(config):
    dna_region = config.dna_region
    config.dna_region = DNARegionType.genic

    fn_pkl = 'dict_cpg_gene.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_cpg_gene = pickle.load(f)
        f.close()
    else:
        dict_cpg_gene = {}

        annotations = config.annotations

        cpg_names = annotations[Annotation.cpg.value]
        gene_names = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        geo = annotations[Annotation.geo.value]
        map_info = annotations[Annotation.map_info.value]
        snp = annotations[Annotation.Probe_SNPs.value]
        snp_10 = annotations[Annotation.Probe_SNPs_10.value]
        cross_reactive = annotations[Annotation.cross_reactive.value]
        for i in range(0, len(cpg_names)):

            curr_cpg = cpg_names[i]
            curr_gene = gene_names[i]
            curr_chr = chr[i]
            curr_geo = geo[i]
            curr_map_info = map_info[i]
            curr_snp = snp[i]
            curr_snp_10 = snp_10[i]
            curr_cross_reactive = cross_reactive[i]

            annotation = {}
            annotation[Annotation.cpg.value] = curr_cpg
            annotation[Annotation.gene.value] = curr_gene
            annotation[Annotation.chr.value] = curr_chr
            annotation[Annotation.geo.value] = curr_geo
            annotation[Annotation.map_info.value] = curr_map_info
            annotation[Annotation.Probe_SNPs.value] = curr_snp
            annotation[Annotation.Probe_SNPs_10.value] = curr_snp_10
            annotation[Annotation.cross_reactive.value] = curr_cross_reactive

            if gene_condition(config, annotation):
                all_genes = list(set(curr_gene.split(';')))
                dict_cpg_gene[curr_cpg] = all_genes

        # cross_reactive strict checking
        if config.cross_reactive is CrossReactiveType.cross_reactive_excluded:
            cpgs_for_del = []
            for curr_cpg in dict_cpg_gene:
                cpg_index = cpg_names.index(curr_cpg)
                curr_cross_reactive = cross_reactive[cpg_index]
                if curr_cross_reactive == 1:
                    cpgs_for_del.append(curr_cpg)
            for curr_cpg in cpgs_for_del:
                del dict_cpg_gene[curr_cpg]

        # snp strict checking
        if config.snp is SNPType.snp_excluded:
            cpgs_for_del = []
            for curr_cpg in dict_cpg_gene:
                cpg_index = cpg_names.index(curr_cpg)
                curr_snp = snp[cpg_index]
                curr_snp_10 = snp_10[cpg_index]
                if curr_snp != '' or curr_snp_10 != '':
                    cpgs_for_del.append(curr_cpg)
            for curr_cpg in cpgs_for_del:
                del dict_cpg_gene[curr_cpg]

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_cpg_gene, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    config.dna_region = dna_region

    return dict_cpg_gene


def get_dict_gene_chr(config):
    dna_region = config.dna_region
    config.dna_region = DNARegionType.genic

    fn_pkl = 'dict_gene_chr.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_gene_chr = pickle.load(f)
        f.close()
    else:
        dict_gene_chr = {}

        annotations = config.annotations

        cpg_names = annotations[Annotation.cpg.value]
        gene_names = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        geo = annotations[Annotation.geo.value]
        map_info = annotations[Annotation.map_info.value]
        snp = annotations[Annotation.Probe_SNPs.value]
        snp1_10 = annotations[Annotation.Probe_SNPs_10.value]
        cross_reactive = annotations[Annotation.cross_reactive.value]
        for i in range(0, len(cpg_names)):

            curr_cpg = cpg_names[i]
            curr_gene = gene_names[i]
            curr_chr = chr[i]
            curr_geo = geo[i]
            curr_map_info = map_info[i]
            curr_snp = snp[i]
            curr_snp_10 = snp1_10[i]
            curr_cross_reactive = cross_reactive[i]

            annotation = {}
            annotation[Annotation.cpg.value] = curr_cpg
            annotation[Annotation.gene.value] = curr_gene
            annotation[Annotation.chr.value] = curr_chr
            annotation[Annotation.geo.value] = curr_geo
            annotation[Annotation.map_info.value] = curr_map_info
            annotation[Annotation.Probe_SNPs.value] = curr_snp
            annotation[Annotation.Probe_SNPs_10.value] = curr_snp_10
            annotation[Annotation.cross_reactive.value] = curr_cross_reactive

            if gene_condition(config, annotation):
                all_genes = list(set(curr_gene.split(';')))
                for g in all_genes:
                    if g in dict_gene_chr:
                        if curr_chr not in dict_gene_chr[g]:
                            dict_gene_chr[g].append(curr_chr)
                    else:
                        dict_gene_chr[g] = [curr_chr]

        # cross_reactive strict checking
        if config.cross_reactive is CrossReactiveType.cross_reactive_excluded:
            genes_for_del = []
            dict_gene_cpgs = get_dict_gene_cpg(config)
            for curr_gene in dict_gene_chr:
                if curr_gene not in dict_gene_cpgs:
                    genes_for_del.append(curr_gene)
            for curr_gene in genes_for_del:
                del dict_gene_chr[curr_gene]

        # cross_reactive strict checking
        if config.snp is SNPType.snp_excluded:
            genes_for_del = []
            dict_gene_cpgs = get_dict_gene_cpg(config)
            for curr_gene in dict_gene_chr:
                if curr_gene not in dict_gene_cpgs:
                    genes_for_del.append(curr_gene)
            for curr_gene in genes_for_del:
                del dict_gene_chr[curr_gene]

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_gene_chr, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    config.dna_region = dna_region

    return dict_gene_chr


def get_dict_gene_cpg(config):
    dna_region = config.dna_region
    config.dna_region = DNARegionType.genic

    fn_pkl = 'dict_gene_cpg.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_gene_cpg = pickle.load(f)
        f.close()
    else:
        dict_gene_cpg = {}

        annotations = config.annotations

        cpg_names = annotations[Annotation.cpg.value]
        gene_names = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        geo = annotations[Annotation.geo.value]
        map_info = annotations[Annotation.map_info.value]
        snp = annotations[Annotation.Probe_SNPs.value]
        snp_10 = annotations[Annotation.Probe_SNPs_10.value]
        cross_reactive = annotations[Annotation.cross_reactive.value]
        for i in range(0, len(cpg_names)):

            curr_cpg = cpg_names[i]
            curr_gene = gene_names[i]
            curr_chr = chr[i]
            curr_geo = geo[i]
            curr_map_info = map_info[i]
            curr_snp = snp[i]
            curr_snp_10 = snp_10[i]
            curr_cross_reactive = cross_reactive[i]

            annotation = {}
            annotation[Annotation.cpg.value] = curr_cpg
            annotation[Annotation.gene.value] = curr_gene
            annotation[Annotation.chr.value] = curr_chr
            annotation[Annotation.geo.value] = curr_geo
            annotation[Annotation.map_info.value] = curr_map_info
            annotation[Annotation.Probe_SNPs.value] = curr_snp
            annotation[Annotation.Probe_SNPs_10.value] = curr_snp_10
            annotation[Annotation.cross_reactive.value] = curr_cross_reactive

            if gene_condition(config, annotation):
                all_genes = list(set(curr_gene.split(';')))
                for g in all_genes:
                    if g in dict_gene_cpg:
                        dict_gene_cpg[g].append(curr_cpg)
                    else:
                        dict_gene_cpg[g] = [curr_cpg]

        # Sorting cpgs in bops by map_info
        for curr_gene, curr_cpgs in dict_gene_cpg.items():
            map_infos = []
            for curr_cpg in curr_cpgs:
                cpg_index = cpg_names.index(curr_cpg)
                curr_map_info = map_info[cpg_index]
                map_infos.append(curr_map_info)
            order = np.argsort(map_infos)
            cpgs_sorted = list(np.array(curr_cpgs)[order])
            dict_gene_cpg[curr_gene] = cpgs_sorted

        # cross_reactive strict checking
        if config.cross_reactive is CrossReactiveType.cross_reactive_excluded:
            genes_for_del = []
            for curr_gene, curr_cpgs in dict_gene_cpg.items():
                for curr_cpg in curr_cpgs:
                    cpg_index = cpg_names.index(curr_cpg)
                    curr_cross_reactive = cross_reactive[cpg_index]
                    if curr_cross_reactive == 1:
                        genes_for_del.append(curr_gene)
                        break
            for curr_gene in genes_for_del:
                del dict_gene_cpg[curr_gene]

        # snp strict checking
        if config.snp is SNPType.snp_excluded:
            genes_for_del = []
            for curr_gene, curr_cpgs in dict_gene_cpg.items():
                for curr_cpg in curr_cpgs:
                    cpg_index = cpg_names.index(curr_cpg)
                    curr_snp = snp[cpg_index]
                    curr_snp_10 = snp_10[cpg_index]
                    if curr_snp != '' or curr_snp_10 != '':
                        genes_for_del.append(curr_gene)
                        break
            for curr_gene in genes_for_del:
                del dict_gene_cpg[curr_gene]

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_gene_cpg, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    config.dna_region = dna_region

    return dict_gene_cpg


def get_dict_cpg_map_info(config):
    dna_region = config.dna_region
    config.dna_region = DNARegionType.genic

    fn_pkl = 'dict_cpg_map_info.pkl'
    fn_pkl = get_data_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_cpg_map_info = pickle.load(f)
        f.close()
    else:
        dict_cpg_map_info = {}

        annotations = config.annotations

        cpg_names = annotations[Annotation.cpg.value]
        gene_names = annotations[Annotation.gene.value]
        chr = annotations[Annotation.chr.value]
        geo = annotations[Annotation.geo.value]
        map_info = annotations[Annotation.map_info.value]
        snp = annotations[Annotation.Probe_SNPs.value]
        snp_10 = annotations[Annotation.Probe_SNPs_10.value]
        cross_reactive = annotations[Annotation.cross_reactive.value]
        for i in range(0, len(cpg_names)):

            curr_cpg = cpg_names[i]
            curr_gene = gene_names[i]
            curr_chr = chr[i]
            curr_geo = geo[i]
            curr_map_info = map_info[i]
            curr_snp = snp[i]
            curr_snp_10 = snp_10[i]
            curr_cross_reactive = cross_reactive[i]

            annotation = {}
            annotation[Annotation.cpg.value] = curr_cpg
            annotation[Annotation.gene.value] = curr_gene
            annotation[Annotation.chr.value] = curr_chr
            annotation[Annotation.geo.value] = curr_geo
            annotation[Annotation.map_info.value] = curr_map_info
            annotation[Annotation.Probe_SNPs.value] = curr_snp
            annotation[Annotation.Probe_SNPs_10.value] = curr_snp_10
            annotation[Annotation.cross_reactive.value] = curr_cross_reactive

            if gene_condition(config, annotation):
                dict_cpg_map_info[curr_cpg] = curr_map_info

        # cross_reactive strict checking
        if config.cross_reactive is CrossReactiveType.cross_reactive_excluded:
            cpgs_for_del = []
            for curr_cpg in dict_cpg_map_info:
                cpg_index = cpg_names.index(curr_cpg)
                curr_cross_reactive = cross_reactive[cpg_index]
                if curr_cross_reactive == 1:
                    cpgs_for_del.append(curr_cpg)
            for curr_cpg in cpgs_for_del:
                del dict_cpg_map_info[curr_cpg]

        # snp strict checking
        if config.snp is SNPType.snp_excluded:
            cpgs_for_del = []
            for curr_cpg in dict_cpg_map_info:
                cpg_index = cpg_names.index(curr_cpg)
                curr_snp = snp[cpg_index]
                curr_snp_10 = snp_10[cpg_index]
                if curr_snp != '' or curr_snp_10 != '':
                    cpgs_for_del.append(curr_cpg)
            for curr_cpg in cpgs_for_del:
                del dict_cpg_map_info[curr_cpg]

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_cpg_map_info, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    config.dna_region = dna_region

    return dict_cpg_map_info
