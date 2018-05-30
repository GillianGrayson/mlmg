import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dicts

type = FSType.local_msi
num_top = 100

suffix = '_islands'

fn = 'table.txt'
full_path = get_full_path(type, fn)
file = open(full_path)
table = file.read().splitlines()

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
genes_pval_mean = []
genes_pval_std = []
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
    genes_pval_mean.append(anova_mean.pvalue)

    anova_std = stats.f_oneway(*std_dict.values())
    genes_pval_std.append(anova_std.pvalue)

order_mean = np.argsort(genes_pval_mean)
pvals_opt_mean = list(np.array(genes_pval_mean)[order_mean])
genes_opt_mean = list(np.array(genes)[order_mean])
info = np.zeros(len(list(genes_opt_mean)), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_mean
info['var2'] = pvals_opt_mean
np.savetxt('anova_genes_mean.txt', info, fmt=fmt)

order_std = np.argsort(genes_pval_std)
pvals_opt_std = list(np.array(genes_pval_std)[order_std])
genes_opt_std = list(np.array(genes)[order_std])
pvals_opt_std_from = list(np.array(genes_pval_std)[order_mean])
genes_opt_std_from = list(np.array(genes)[order_mean])
info = np.zeros(len(list(genes_opt_std)), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_std
info['var2'] = pvals_opt_std
np.savetxt('anova_genes_std.txt', info, fmt=fmt)
info = np.zeros(len(list(pvals_opt_std_from)), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_std_from
info['var2'] = pvals_opt_std_from
np.savetxt('anova_genes_std_from.txt', info, fmt=fmt)

genes_match = []
for gene in genes_opt_mean[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('anova_match_mean.txt', genes_match, fmt="%s")
print('top: ' + str(len(genes_match)))

genes_match = []
for gene in genes_opt_std[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('anova_match_std.txt', genes_match, fmt="%s")
print('top: ' + str(len(genes_match)))
