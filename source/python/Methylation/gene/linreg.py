import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dicts
from gen_files.geo import *

fs_type = FSType.local_big
db_type = DataBaseType.GSE52588
geo_type = GeoType.islands_shores

fn = 'attribute.txt'
ages = []
full_path = get_full_path(fs_type, db_type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

sort_type = 1 # 0 - pval, 1 - rval

mean_names = []
mean = []
mean_p_values = []
mean_r_values = []
mean_slopes = []
mean_intercepts = []

std_names = []
std = []
std_p_values = []
std_r_values = []
std_slopes = []
std_intercepts = []

mean_der_names = []
mean_der = []
mean_der_p_values = []
mean_der_r_values = []
mean_der_slopes = []
mean_der_intercepts = []

fn = 'gene_mean' + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    mean_names.append(gene)
    mean.append(vals)

fn = 'gene_std' + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    std_names.append(gene)
    std.append(vals)

fn = 'gene_mean_der' + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    mean_der_names.append(gene)
    mean_der.append(vals)

for id in range(0, len(mean_names)):

    means = mean[id]
    stds = std[id]
    mean_ders = mean_der[id]

    gene_name_mean = mean_names[id]
    gene_name_std = std_names[id]
    gene_name_mean_der = mean_der_names[id]

    if gene_name_mean != gene_name_std:
        print('error')

    slope, intercept, r_value, p_value, std_err = stats.linregress(means, ages)
    mean_r_values.append(r_value)
    mean_p_values.append(p_value)
    mean_slopes.append(slope)
    mean_intercepts.append(intercept)

    slope, intercept, r_value, p_value, std_err = stats.linregress(stds, ages)
    std_r_values.append(r_value)
    std_p_values.append(p_value)
    std_slopes.append(slope)
    std_intercepts.append(intercept)

    slope, intercept, r_value, p_value, std_err = stats.linregress(mean_ders, ages)
    mean_der_r_values.append(r_value)
    mean_der_p_values.append(p_value)
    mean_der_slopes.append(slope)
    mean_der_intercepts.append(intercept)

order_mean= []
if sort_type == 0:
    order_mean = np.argsort(mean_p_values)
else:
    order_mean = np.argsort(list(map(abs, mean_r_values)))[::-1]
p_values_opt_mean = list(np.array(mean_p_values)[order_mean])
r_values_opt_mean = list(np.array(mean_r_values)[order_mean])
slopes_opt_mean = list(np.array(mean_slopes)[order_mean])
intercepts_opt_mean = list(np.array(mean_intercepts)[order_mean])
genes_opt_mean = list(np.array(mean_names)[order_mean])
info = np.zeros(len(mean_names), dtype=[('var1', 'U50'), ('var2', float), ('var3', float), ('var4', float), ('var5', float)])
fmt = "%s %0.18e %0.18e %0.18e %0.18e"
info['var1'] = genes_opt_mean
info['var2'] = p_values_opt_mean
info['var3'] = r_values_opt_mean
info['var4'] = slopes_opt_mean
info['var5'] = intercepts_opt_mean
np.savetxt('linreg_genes_mean' + geo_type.value + '.txt', info, fmt=fmt)

order_std= []
if sort_type == 0:
    order_std = np.argsort(std_p_values)
else:
    order_std = np.argsort(list(map(abs, std_r_values)))[::-1]
p_values_opt_std = list(np.array(std_p_values)[order_std])
r_values_opt_std = list(np.array(std_r_values)[order_std])
slopes_opt_std = list(np.array(std_slopes)[order_std])
intercepts_opt_std = list(np.array(std_intercepts)[order_std])
genes_opt_std = list(np.array(std_names)[order_std])
info = np.zeros(len(std_names), dtype=[('var1', 'U50'), ('var2', float), ('var3', float), ('var4', float), ('var5', float)])
fmt = "%s %0.18e %0.18e %0.18e %0.18e"
info['var1'] = genes_opt_std
info['var2'] = p_values_opt_std
info['var3'] = r_values_opt_std
info['var4'] = slopes_opt_std
info['var5'] = intercepts_opt_std
np.savetxt('linreg_genes_std' + geo_type.value + '.txt', info, fmt=fmt)

order_mean_der= []
if sort_type == 0:
    order_mean_der = np.argsort(mean_der_p_values)
else:
    order_mean_der = np.argsort(list(map(abs, mean_der_r_values)))[::-1]
p_values_opt_mean_der = list(np.array(mean_der_p_values)[order_mean_der])
r_values_opt_mean_der = list(np.array(mean_der_r_values)[order_mean_der])
slopes_opt_mean_der = list(np.array(mean_der_slopes)[order_mean_der])
intercepts_opt_mean_der = list(np.array(mean_der_intercepts)[order_mean_der])
genes_opt_mean_der = list(np.array(mean_der_names)[order_mean_der])
info = np.zeros(len(mean_der_names), dtype=[('var1', 'U50'), ('var2', float), ('var3', float), ('var4', float), ('var5', float)])
fmt = "%s %0.18e %0.18e %0.18e %0.18e"
info['var1'] = genes_opt_mean_der
info['var2'] = p_values_opt_mean_der
info['var3'] = r_values_opt_mean_der
info['var4'] = slopes_opt_mean_der
info['var5'] = intercepts_opt_mean_der
np.savetxt('linreg_genes_mean_der' + geo_type.value + '.txt', info, fmt=fmt)

if db_type is DataBaseType.GSE40279:

    num_top = 100

    fn = 'table.txt'
    full_path = get_full_path(fs_type, db_type, fn)
    file = open(full_path)
    table = file.read().splitlines()

    genes_match = []
    for gene in genes_opt_mean[0:num_top]:
        if gene in table:
            genes_match.append(gene)
    np.savetxt('linreg_match_mean.txt', genes_match, fmt="%s")
    print('top: ' + str(len(genes_match)))


    genes_match = []
    for gene in genes_opt_std[0:num_top]:
        if gene in table:
            genes_match.append(gene)
    np.savetxt('linreg_match_std.txt', genes_match, fmt="%s")
    print('top: ' + str(len(genes_match)))