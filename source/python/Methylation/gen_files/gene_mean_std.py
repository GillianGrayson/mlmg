import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import *
from gen_files.geo import *
from config import *

print_rate = 10000

fs_type = FSType.local_big
db_type = DataBaseType.GSE52588
geo_type = GeoType.islands_shores
config = Config(fs_type, db_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type)

dict_cpg_gene = get_dict_cpg_gene(fs_type, db_type, geo_type)

gene_raw_dict = config.get_raw_dict(dict_cpg_gene)

gene_mean_dict = {}
gene_std_dict = {}
gene_mean_str_list = []
gene_std_str_list = []
num_genes = 0
for gene in gene_raw_dict:

    mean_list = []
    std_list = []
    for curr_list in gene_raw_dict[gene]:
        mean_list.append(np.mean(curr_list))
        std_list.append(np.std(curr_list))

    gene_mean_dict[gene] = mean_list
    gene_std_dict[gene] = std_list

    curr_mean_str = gene
    curr_std_str = gene
    for id in range(0, len(mean_list)):
        curr_mean_str += (' ' + str(format(mean_list[id], '0.8e')))
        curr_std_str += (' ' + str(format(std_list[id], '0.8e')))

    gene_mean_str_list.append(curr_mean_str)
    gene_std_str_list.append(curr_std_str)

    num_genes += 1
    if num_genes % print_rate == 0:
        print('num_genes: ' + str(num_genes))

np.savetxt('gene_mean' + geo_type.value + '.txt', gene_mean_str_list, fmt="%s")
np.savetxt('gene_std' + geo_type.value + '.txt', gene_std_str_list, fmt="%s")


