import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import *

type = FSType.local_msi
print_rate = 1000

suffix = '_islands_shores'

dict_cpg_gene = {}
if suffix == '':
    dict_cpg_gene = get_dicts(type)
elif suffix == '_shores':
    dict_cpg_gene = get_dict_cpg_gene_shore(type)
elif suffix == '_islands':
    dict_cpg_gene = get_dict_cpg_gene_island(type)
elif suffix == '_islands_shores':
    dict_cpg_gene = get_dict_cpg_gene_island_shore(type)


fn = 'ages.txt'
ages = []
full_path = get_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

fn = 'GSE40279_average_beta.txt'
full_path = get_path(type, fn)

num_lines = 0

cpgs = []
rhos = []
pvals = []

f = open(full_path)
first_line = f.readline()
col_names = first_line.split('\t')

for line in f:

    col_vals = line.split('\t')
    CpG = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    if CpG in dict_cpg_gene:

        rho, pval = stats.spearmanr(ages, vals)

        genes = dict_cpg_gene.get(CpG)

        if genes is not None:

            cpgs.append(CpG)
            rhos.append(rho)
            pvals.append(pval)

            num_lines += 1

            if num_lines % print_rate == 0:
                print('num_lines: ' + str(num_lines))

order = np.argsort(list(map(abs, rhos)))[::-1]
rhos_opt = list(np.array(rhos)[order])
cpgs_opt = list(np.array(cpgs)[order])
pvals_opt = list(np.array(pvals)[order])

genes_spec = []
rhos_spec = []
for id in range(0, len(cpgs_opt)):
    cpg = cpgs_opt[id]
    rho = rhos_opt[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        if gene not in genes_spec:
            genes_spec.append(gene)
            rhos_spec.append(rho)

info = np.zeros(len(genes_spec), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(genes_spec)
info['var2'] = list(rhos_spec)
np.savetxt('spearman_spec.txt', info, fmt=fmt)

genes_names_mean = []
genes_names_std = []
genes_mean = []
genes_std = []
genes_mean_rho = {}
genes_std_rho = {}

fn = 'gene_mean' + suffix + '.txt'
full_path = get_path(type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes_names_mean.append(gene)
    genes_mean.append(vals)

fn = 'gene_std' + suffix + '.txt'
full_path = get_path(type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes_names_std.append(gene)
    genes_std.append(vals)

for id in range(0, len(genes_names_mean)):

    means = genes_mean[id]
    stds = genes_std[id]

    gene_name_mean = genes_names_mean[id]
    gene_name_std = genes_names_std[id]

    if gene_name_mean != gene_name_std:
        print('error')

    rho_means, pval_means = stats.spearmanr(ages, means)
    rho_stds, pval_stds = stats.spearmanr(ages, stds)
    if math.isnan(rho_stds):
        rho_stds = 0.0

    genes_mean_rho[gene_name_mean] = rho_means
    genes_std_rho[gene_name_std] = rho_stds

rhos_spec_from_mean = []
rhos_spec_from_std = []
for gene in genes_spec:
    rhos_spec_from_mean.append(genes_mean_rho[gene])
    rhos_spec_from_std.append(genes_std_rho[gene])

info = np.zeros(len(genes_spec), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(genes_spec)
info['var2'] = list(rhos_spec_from_mean)
np.savetxt('spearman_spec_from_mean.txt', info, fmt=fmt)

info = np.zeros(len(genes_spec), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(genes_spec)
info['var2'] = list(rhos_spec_from_std)
np.savetxt('spearman_spec_from_std.txt', info, fmt=fmt)