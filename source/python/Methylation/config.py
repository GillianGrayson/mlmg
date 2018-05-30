from Infrastructure.file_system import *
from gen_files.geo import *
from dicts import *
import numpy as np

class Config:

    def __init__(self, fs_type, db_type):
        self.fs_type = fs_type
        self.db_type = db_type
        self.print_rate = 10000

    def get_raw_dict(self, fs_type, db_type, geo_type):
        raise NotImplementedError("Subclass must implement abstract method")


class ConfigGSE40279(Config):

    def __init__(self, fs_type, db_type):
        super().__init__(fs_type, db_type)
        self.num_skip_lines = 1
        self.attribute_fn = 'attribute.txt'

    def get_raw_dict(self, fs_type, db_type, geo_type):

        dict_cpg_gene, dict_cpg_map = get_dicts(fs_type, db_type, geo_type)

        fn = self.attribute_fn
        attribute = []
        full_path = get_full_path(self.fs_type, self.db_type, fn)
        with open(full_path) as f:
            for line in f:
                attribute.append(int(line))

        fn = self.db_type.value + '_average_beta.txt'
        full_path = get_full_path(self.fs_type, self.db_type, fn)
        f = open(full_path)
        for skip_id in range(0, self.num_skip_lines):
            skip_line = f.readline()

        num_lines = 0
        gene_raw_dict = {}
        map_dict = {}
        for line in f:

            col_vals = line.split('\t')
            CpG = col_vals[0]
            vals = list(map(float, col_vals[1::]))

            genes = dict_cpg_gene.get(CpG)
            map_info = dict_cpg_map.get(CpG)

            if genes is not None:
                for gene in genes:
                    if gene in gene_raw_dict:
                        for list_id in range(0, len(attribute)):
                            gene_raw_dict[gene][list_id].append(vals[list_id])
                        map_dict[gene].append(int(map_info))
                    else:
                        gene_raw_dict[gene] = []
                        for list_id in range(0, len(attribute)):
                            gene_raw_dict[gene].append([vals[list_id]])
                        map_dict[gene] = []
                        map_dict[gene].append(int(map_info))

            num_lines += 1
            if num_lines % self.print_rate == 0:
                print('num_lines: ' + str(num_lines))

        for gene in gene_raw_dict:
            raw = gene_raw_dict[gene]
            map_info = map_dict[gene]
            order = np.argsort(map_info)
            gene_raw_dict[gene] = []
            for record in raw:
                sorted_record = list(np.array(record)[order])
                gene_raw_dict[gene].append(sorted_record)


        return gene_raw_dict



class ConfigGSE52588(Config):
    def __init__(self, fs_type, db_type):
        super().__init__(fs_type, db_type)
        self.num_skip_lines = 87
        self.attribute_fn = 'attribute.txt'

    def get_raw_dict(self, fs_type, db_type, geo_type):

        dict_cpg_gene, dict_cpg_map = get_dicts(fs_type, db_type, geo_type)

        num_skip_lines_raw = 1
        pval_lim = 0.05
        pval_part = 0.75

        cpg_non_inc = []
        fn = self.db_type.value + '_raw_data.txt'
        full_path = get_full_path(self.fs_type, self.db_type, fn)
        f = open(full_path)
        for skip_id in range(0, num_skip_lines_raw):
            skip_line = f.readline()

        num_lines = 0
        for line in f:
            col_vals = line.split('\t')
            CpG = col_vals[0].replace('"', '')
            vals = list(map(float, col_vals[1::]))
            pvals = vals[2::3]
            num_big_pvals = 0

            for pval in pvals:
                if pval > pval_lim:
                    num_big_pvals += 1

            if float(num_big_pvals) / float(len(pvals)) > pval_part:
                cpg_non_inc.append(CpG)

            num_lines += 1
            if num_lines % self.print_rate == 0:
                print('num_lines: ' + str(num_lines))
                print('num_cpg_non_inc: ' + str(len(cpg_non_inc)))

        fn = self.attribute_fn
        attribute = []
        full_path = get_full_path(self.fs_type, self.db_type, fn)
        with open(full_path) as f:
            for line in f:
                attribute.append(int(line))

        fn = self.db_type.value + '_average_beta.txt'
        full_path = get_full_path(self.fs_type, self.db_type, fn)
        f = open(full_path)
        for skip_id in range(0, self.num_skip_lines):
            skip_line = f.readline()

        num_miss = 0
        num_lines = 0
        gene_raw_dict = {}
        map_dict = {}
        for line in f:

            col_vals = line.split('\t')

            is_none = False
            for val_id in range(0, len(col_vals)):
                col_vals[val_id] = col_vals[val_id].replace('"', '').rstrip()

            if 'NULL' in col_vals:
                is_none = True
                num_miss += 1

            if not is_none:

                CpG = col_vals[0].replace('"', '')

                if CpG not in cpg_non_inc:

                    vals = list(map(float, col_vals[1::]))

                    genes = dict_cpg_gene.get(CpG)
                    map_info = dict_cpg_map.get(CpG)

                    if genes is not None:
                        for gene in genes:
                            if gene in gene_raw_dict:
                                for list_id in range(0, len(attribute)):
                                    gene_raw_dict[gene][list_id].append(vals[list_id])
                                map_dict[gene].append(int(map_info))
                            else:
                                gene_raw_dict[gene] = []
                                for list_id in range(0, len(attribute)):
                                    gene_raw_dict[gene].append([vals[list_id]])
                                map_dict[gene] = []
                                map_dict[gene].append(int(map_info))

                    num_lines += 1
                    if num_lines % self.print_rate == 0:
                        print('num_lines: ' + str(num_lines))

        print('num_miss: ' + str(num_miss))
        print('num_final: ' + str(num_lines))

        for gene in gene_raw_dict:
            raw = gene_raw_dict[gene]
            map_info = map_dict[gene]
            order = np.argsort(map_info)
            gene_raw_dict[gene] = []
            for record in raw:
                sorted_record = list(np.array(record)[order])
                gene_raw_dict[gene].append(sorted_record)

        return gene_raw_dict