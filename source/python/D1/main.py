import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from cpg_gene_dict import get_dict_cpg_gene

type = FSType.local
pval_lim = 1.0e-10
print_rate = 1000
num_pval_genes = 1000

dict_cpg_gene = get_dict_cpg_gene(type)

fn = 'ages.txt'
ages = []
full_path = get_full_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

shift_age = 10
min_age = min(ages)
max_age = max(ages)

min_age = int(min_age / shift_age) * shift_age
max_age = (int(max_age / shift_age) + 1) * shift_age

fn = 'GSE40279_average_beta.txt'
full_path = get_full_path(type, fn)

age_dict = {}
for age_id in range(0, len(ages)):
    age = ages[age_id]
    key = int((age - min_age) / shift_age)
    if key in age_dict:
        age_dict[key].append(age_id)
    else:
        age_dict[key] = [age_id]

num_lines = 0
approved_CpGs = []
approved_pvals = []
approved_genes = []
genes_rate_dict = {}
f = open(full_path)
first_line = f.readline()
col_names = first_line.split('\t')

for line in f:

    col_vals = line.split('\t')
    CpG = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    curr_beta_dict = {}
    for key_age in age_dict:
        curr_beta_dict[key_age] = list(np.asarray(vals)[age_dict[key_age]])

    anova_res = stats.f_oneway(*curr_beta_dict.values())

    if(anova_res.pvalue < pval_lim):

        approved_CpGs.append(CpG)
        approved_pvals.append(anova_res.pvalue)
        genes = dict_cpg_gene.get(CpG)

        if genes is None:

            approved_genes.append('n')

            if 'n' in genes_rate_dict:
                genes_rate_dict['n'] += 1
            else:
                genes_rate_dict['n'] = 1


        else:

            for gene in genes:
                if gene in genes_rate_dict:
                    genes_rate_dict[gene] += 1
                else:
                    genes_rate_dict[gene] = 1

            approved_genes.append(';'.join(genes))

    num_lines += 1

    if num_lines % print_rate == 0:
        print('num_lines: ' + str(num_lines))
        print('num_CpGs: ' + str(len(approved_CpGs)))

info = np.zeros(len(approved_CpGs), dtype=[('var1', 'U50'), ('var2', float), ('var3', 'U50')])
fmt = "%s %18e %s"
info['var1'] = approved_CpGs
info['var2'] = approved_pvals
info['var3'] = approved_genes
np.savetxt('approved_CpGs.txt', info, fmt=fmt)

info = np.zeros(len(list(genes_rate_dict.keys())), dtype=[('var1', 'U50'), ('var2', int)])
fmt = "%s %d"
info['var1'] = list(genes_rate_dict.keys())
info['var2'] = list(genes_rate_dict.values())
np.savetxt('genes_rate.txt', info, fmt=fmt)

order = np.argsort(approved_pvals)
min_pvals = sorted(approved_pvals)[0:num_pval_genes]
min_genes = list(np.array(approved_genes)[order[0:num_pval_genes]])
min_genes_dict = {}
for i in range(0, len(min_pvals)):
    curr_genes = min_genes[i]
    curr_pval = min_pvals[i]
    curr_min_genes = list(set(curr_genes.split(';')))
    for gene in curr_min_genes:
        if gene in min_genes_dict:
            min_genes_dict[gene] += 1
        else:
            min_genes_dict[gene] = 1

info = np.zeros(len(list(min_genes_dict.keys())), dtype=[('var1', 'U50'), ('var2', int)])
fmt = "%s %d"
info['var1'] = list(min_genes_dict.keys())
info['var2'] = list(min_genes_dict.values())
np.savetxt('pvals_genes.txt', info, fmt=fmt)
