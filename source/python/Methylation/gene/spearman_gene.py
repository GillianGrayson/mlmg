import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from gen_files.geo import *
from dicts import get_dicts

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.islands_shores
num_top = 100

fn = 'table.txt'
full_path = get_full_path(fs_type, db_type, fn)
file = open(full_path)
table = file.read().splitlines()

fn = 'attribute.txt'
ages = []
full_path = get_full_path(fs_type, db_type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

num_genes = 0
genes_names_mean = []
genes_names_std = []
genes_mean = []
genes_std = []
genes_mean_rho = []
genes_std_rho = []

fn = 'gene_mean' + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes_names_mean.append(gene)
    genes_mean.append(vals)

fn = 'gene_std' + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
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

    genes_mean_rho.append(rho_means)
    genes_std_rho.append(rho_stds)

order_mean = np.argsort(list(map(abs, genes_mean_rho)))[::-1]
rho_opt_mean = list(np.array(genes_mean_rho)[order_mean])
genes_opt_mean = list(np.array(genes_names_mean)[order_mean])
info = np.zeros(len(genes_names_mean), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_mean
info['var2'] = rho_opt_mean
np.savetxt('spearman_genes_mean.txt', info, fmt=fmt)

order_std = np.argsort(list(map(abs, genes_std_rho)))[::-1]
rho_opt_std = list(np.array(genes_std_rho)[order_std])
genes_opt_std = list(np.array(genes_names_std)[order_std])
rho_opt_std_from = list(np.array(genes_std_rho)[order_mean])
genes_opt_std_from = list(np.array(genes_names_std)[order_mean])
info = np.zeros(len(genes_names_std), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_std
info['var2'] = rho_opt_std
np.savetxt('spearman_genes_std.txt', info, fmt=fmt)
info = np.zeros(len(genes_names_std), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = genes_opt_std_from
info['var2'] = rho_opt_std_from
np.savetxt('spearman_genes_std_from.txt', info, fmt=fmt)

genes_match = []
for gene in genes_opt_mean[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('spearman_match_mean.txt', genes_match, fmt="%s")
print('top: ' + str(len(genes_match)))


genes_match = []
for gene in genes_opt_std[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('spearman_match_std.txt', genes_match, fmt="%s")
print('top: ' + str(len(genes_match)))