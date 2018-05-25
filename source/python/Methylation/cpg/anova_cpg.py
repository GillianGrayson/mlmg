import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dict_cpg_gene

type = FSType.local_big
print_rate = 1000
num_opt = 1000
num_top = 100

suffix = '_ws'

fn = 'table.txt'
full_path = get_full_path(type, fn)
file = open(full_path)
table = file.read().splitlines()

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
cpgs = []
pvals = []
f = open(full_path)
first_line = f.readline()
col_names = first_line.split('\t')

for line in f:

    col_vals = line.split('\t')
    CpG = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes = dict_cpg_gene.get(CpG)

    curr_beta_dict = {}
    for key_age in age_dict:
        curr_beta_dict[key_age] = list(np.asarray(vals)[age_dict[key_age]])

    anova_res = stats.f_oneway(*curr_beta_dict.values())

    if genes is not None:

        cpgs.append(CpG)
        pvals.append(anova_res.pvalue)

    num_lines += 1

    if num_lines % print_rate == 0:
        print('num_lines: ' + str(num_lines))

order = np.argsort(pvals)
cpgs_opt = list(np.array(cpgs)[order[0:num_opt]])
pvals_opt = list(np.array(pvals)[order[0:num_opt]])
genes_opt = []
for cpg in cpgs_opt:
    genes = dict_cpg_gene.get(cpg)
    sum = ''
    for gene in genes:
        sum += str(gene)
    genes_opt.append(sum)


info = np.zeros(len(cpgs_opt), dtype=[('var1', 'U50'), ('var2', 'U50'), ('var3', float)])
fmt = "%s %s %0.18e"
info['var1'] = list(cpgs_opt)
info['var2'] = list(genes_opt)
info['var3'] = list(pvals_opt)
np.savetxt('anova_full.txt', info, fmt=fmt)

genes_set = []
pval_set = []
for gene_id in range(0, len(genes_opt)):
    gene = genes_opt[gene_id]
    pval = pvals_opt[gene_id]
    if gene not in genes_set:
        genes_set.append(gene)
        pval_set.append(pval)
np.savetxt('anova_genes.txt', genes_set[0:num_top], fmt="%s")
np.savetxt('anova_pvals.txt', pval_set[0:num_top], fmt="%0.18e")

genes_match = []
for gene in genes_set[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('anova_match.txt', genes_match, fmt="%s")

print('top: ' + str(len(genes_match)))
