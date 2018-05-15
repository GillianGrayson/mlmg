import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dict_cpg_gene

type = FSType.local_msi
print_rate = 3000

dict_cpg_gene = get_dict_cpg_gene(type)

fn = 'ages.txt'
ages = []
full_path = get_full_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

fn = 'GSE40279_average_beta.txt'
full_path = get_full_path(type, fn)
f = open(full_path)
first_line = f.readline()
col_names = first_line.split('\t')

num_lines = 0
gene_raw_dict = {}
for line in f:

    col_vals = line.split('\t')
    CpG = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    genes = dict_cpg_gene.get(CpG)

    if genes is not None:
        for gene in genes:
            if gene in gene_raw_dict:
                for list_id in range(0, len(ages)):
                    gene_raw_dict[gene][list_id].append(vals[list_id])
            else:
                gene_raw_dict[gene] = []
                for list_id in range(0, len(ages)):
                    gene_raw_dict[gene].append([vals[list_id]])

    num_lines += 1
    if num_lines % print_rate == 0:
        print('num_lines: ' + str(num_lines))
    if num_lines == 1000:
        break

gene_mean_dict = {}
gene_std_dict = {}
for gene in gene_raw_dict:
    mean_list = []
    std_list = []
    for curr_list in gene_raw_dict[gene]:
        mean_list.append(np.mean(curr_list))
        std_list.append(np.std(curr_list))
    gene_mean_dict[gene] = mean_list
    gene_std_dict[gene] = std_list



ololo = 5



