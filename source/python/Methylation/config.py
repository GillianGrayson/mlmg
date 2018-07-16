import numpy as np
from annotation.regular import *
from infrastructure.load import *
from infrastructure.file_system import *

def get_annotations(config):
    fn = 'annotations.txt'
    fn = get_path(config.fs_type, config.db_type, fn)
    f = open(fn)

    key_line = f.readline()
    keys = key_line.split('\t')
    annotations = {}
    for key_id in range(0, len(keys)):
        keys[key_id] = keys[key_id].rstrip()
        annotations[keys[key_id]] = []

    for line in f:
        vals = line.split('\t')
        for key_id in range(0, len(keys)):
            annotations[keys[key_id]].append(vals[key_id].rstrip())

    return annotations

class Config:

    def __init__(self, fs_type, db_type, geo_type=GeoType.any, dna_region = DNARegion.genic, class_type=ClassType.any, clustering_type=ClusteringType.k_means):
        self.fs_type = fs_type
        self.db_type = db_type
        self.geo_type = geo_type
        self.class_type = class_type
        self.dna_region = dna_region
        self.clustering_type = clustering_type
        self.print_rate = 10000

    def get_raw_dict(self):
        raise NotImplementedError("Subclass must implement abstract method")


class ConfigGSE40279(Config):

    def __init__(self, fs_type, db_type, geo_type=GeoType.any, dna_region = DNARegion.genic, class_type=ClassType.any, clustering_type=ClusteringType.k_means):
        super().__init__(fs_type, db_type, geo_type=geo_type, dna_region=dna_region, class_type=class_type, clustering_type=clustering_type)
        self.num_skip_lines = 1
        self.attribute_fn = 'attribute.txt'
        self.miss_tag = 'NULL'
        self.annotations = get_annotations(self)

    def get_raw_dict(self):
        dict_cpg_gene = get_dict_cpg_gene(self)
        dict_cpg_map = get_dict_cpg_map_info(self)
        attributes = get_attributes(self)
        cpgs, vals = get_cpg_data(self)

        gene_raw_dict = {}
        map_dict = {}
        for id in range(0, len(cpgs)):

            curr_cpg = cpgs[id]
            curr_vals = vals[id]

            genes = dict_cpg_gene.get(curr_cpg)
            map_info = dict_cpg_map.get(curr_cpg)

            if genes is not None:
                for gene in genes:
                    if gene in gene_raw_dict:
                        for list_id in range(0, len(attributes)):
                            gene_raw_dict[gene][list_id].append(curr_vals[list_id])
                        map_dict[gene].append(int(map_info))
                    else:
                        gene_raw_dict[gene] = []
                        for list_id in range(0, len(attributes)):
                            gene_raw_dict[gene].append([curr_vals[list_id]])
                        map_dict[gene] = []
                        map_dict[gene].append(int(map_info))

        for gene in gene_raw_dict:
            raw = gene_raw_dict[gene]
            map_info = map_dict[gene]
            order = np.argsort(map_info)
            gene_raw_dict[gene] = []
            for record in raw:
                sorted_record = list(np.array(record)[order])
                gene_raw_dict[gene].append(sorted_record)

        return gene_raw_dict

    def line_proc(self, line):
        line_list = line.split('\t')
        return line_list


class ConfigGSE52588(Config):
    def __init__(self, fs_type, db_type, geo_type=GeoType.any, dna_region = DNARegion.genic, class_type=ClassType.any, clustering_type=ClusteringType.k_means):
        super().__init__(fs_type, db_type, geo_type=geo_type, dna_region=dna_region, class_type=class_type, clustering_type=clustering_type)
        self.num_skip_lines = 87
        self.attribute_fn = 'attribute.txt'
        self.miss_tag = 'NULL'
        self.annotations = get_annotations(self)

    def get_raw_dict(self):
        dict_cpg_gene = get_dict_cpg_gene(self)
        dict_cpg_map = get_dict_cpg_map_info(self)

        pval_lim = 0.05
        pval_part = 0.75

        cpgs, pvals = get_cpg_pval_data(self)
        cpg_non_inc = []
        for id in range(0, len(cpgs)):
            curr_cpg = cpgs[id]
            curr_pvals = pvals[id]

            num_big_pvals = 0
            for pval in curr_pvals:
                if pval > pval_lim:
                    num_big_pvals += 1

            if float(num_big_pvals) / float(len(curr_pvals)) > pval_part:
                cpg_non_inc.append(curr_cpg)

        attributes = get_attributes(self)
        cpgs, vals = get_cpg_data(self)

        gene_raw_dict = {}
        map_dict = {}
        for id in range(0, len(cpgs)):

            curr_cpg = cpgs[id]
            curr_vals = vals[id]

            if curr_cpg not in cpg_non_inc:

                genes = dict_cpg_gene.get(curr_cpg)
                map_info = dict_cpg_map.get(curr_cpg)

                if genes is not None:
                    for gene in genes:
                        if gene in gene_raw_dict:
                            for list_id in range(0, len(attributes)):
                                gene_raw_dict[gene][list_id].append(curr_vals[list_id])
                            map_dict[gene].append(int(map_info))
                        else:
                            gene_raw_dict[gene] = []
                            for list_id in range(0, len(attributes)):
                                gene_raw_dict[gene].append([curr_vals[list_id]])
                            map_dict[gene] = []
                            map_dict[gene].append(int(map_info))

        for gene in gene_raw_dict:
            raw = gene_raw_dict[gene]
            map_info = map_dict[gene]
            order = np.argsort(map_info)
            gene_raw_dict[gene] = []
            for record in raw:
                sorted_record = list(np.array(record)[order])
                gene_raw_dict[gene].append(sorted_record)

        return gene_raw_dict

    def line_proc(self, line):
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()

        return line_list
