import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dict_cpg_gene

type = FSType.local_big
pval_lim = 1.0e-10
print_rate = 1000
num_pval_genes = 1000

suffix = ''

fn = 'ages.txt'
ages = []
full_path = get_full_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

shift_age = 5
min_age = min(ages)
max_age = max(ages)

min_age = int(min_age / shift_age) * shift_age
max_age = (int(max_age / shift_age) + 1) * shift_age
age_dict = {}
for age_id in range(0, len(ages)):
    age = ages[age_id]
    key = int((age - min_age) / shift_age)
    if key in age_dict:
        age_dict[key].append(age_id)
    else:
        age_dict[key] = [age_id]

num_genes = 0
genes = []
genes_mean = []
genes_std = []
genes_mean_pval = []
genes_std_pval = []
genes_mean_std_pval = []

fn = 'gene_mean' + suffix + '.txt'
full_path = get_full_path(type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes.append(gene)
    genes_mean.append(vals)

fn = 'gene_std' + suffix + '.txt'
full_path = get_full_path(type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    vals = list(map(float, col_vals[1::]))
    genes_std.append(vals)

for id in range(0, len(genes)):

    means = genes_mean[id]
    stds = genes_std[id]

    mean_dict = {}
    for key_age in age_dict:
        mean_dict[key_age] = list(np.asarray(means)[age_dict[key_age]])

    std_dict = {}
    for key_age in age_dict:
        std_dict[key_age] = list(np.asarray(stds)[age_dict[key_age]])

    anova_mean = stats.f_oneway(*mean_dict.values())
    genes_mean_pval.append(anova_mean.pvalue)

    anova_std = stats.f_oneway(*std_dict.values())
    genes_std_pval.append(anova_std.pvalue)

order_mean = np.argsort(genes_mean_pval)
min_pvals_mean = sorted(genes_mean_pval)[0:len(genes_mean_pval)]
min_genes_mean = list(np.array(genes)[order_mean[0:len(genes_mean_pval)]])

info = np.zeros(len(list(min_genes_mean)), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = min_genes_mean
info['var2'] = min_pvals_mean
np.savetxt('pvals_mean_genes.txt', info, fmt=fmt)

order_std = np.argsort(genes_std_pval)
min_pvals_std = sorted(genes_std_pval)[0:len(genes_std_pval)]
min_genes_std = list(np.array(genes)[order_mean[0:len(genes_std_pval)]])

info = np.zeros(len(list(min_genes_std)), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = min_genes_std
info['var2'] = min_pvals_std
np.savetxt('pvals_std_genes.txt', info, fmt=fmt)
