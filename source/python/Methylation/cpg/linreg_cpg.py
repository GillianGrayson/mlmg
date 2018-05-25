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

fn = 'table.txt'
full_path = get_full_path(type, fn)
file = open(full_path)
table = file.read().splitlines()

with open(full_path) as f:
    for line in f:
        table.append(str(line))

suffix = ''

dict_cpg_gene = get_dict_cpg_gene(type)

fn = 'ages.txt'
ages = []
full_path = get_full_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

fn = 'GSE40279_average_beta.txt'
full_path = get_full_path(type, fn)

num_lines = 0

cpgs = []
slopes = []
intercepts = []
r_values = []
p_values = []
std_errs = []

f = open(full_path)
first_line = f.readline()
col_names = first_line.split('\t')

for line in f:

    col_vals = line.split('\t')
    CpG = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    if CpG in dict_cpg_gene:

        slope, intercept, r_value, p_value, std_err = stats.linregress(vals, ages)

        genes = dict_cpg_gene.get(CpG)

        if genes is not None:

            cpgs.append(CpG)
            slopes.append(slope)
            intercepts.append(intercept)
            r_values.append(r_value)
            p_values.append(p_value)
            std_errs.append(std_err)

            num_lines += 1

            if num_lines % print_rate == 0:
                print('num_lines: ' + str(num_lines))


order = np.argsort(p_values)

cpgs_opt = list(np.array(cpgs)[order[0:num_opt]])
p_values_opt = list(np.array(p_values)[order[0:num_opt]])
slopes_opt = list(np.array(slopes)[order[0:num_opt]])
intercepts_opt = list(np.array(intercepts)[order[0:num_opt]])
r_values_opt = list(np.array(r_values)[order[0:num_opt]])
std_errs_opt = list(np.array(std_errs)[order[0:num_opt]])

genes_opt = []
for cpg in cpgs_opt:
    genes = dict_cpg_gene.get(cpg)
    sum = ''
    for gene in genes:
        sum += str(gene)
    genes_opt.append(sum)


info = np.zeros(len(cpgs_opt), dtype=[('var1', 'U50'), ('var2', 'U50'), ('var3', float), ('var4', float), ('var5', float), ('var6', float), ('var7', float)])
fmt = "%s %s %0.18e %0.18e %0.18e %0.18e %0.18e"
info['var1'] = list(cpgs_opt)
info['var2'] = list(genes_opt)
info['var3'] = list(p_values_opt)
info['var4'] = list(slopes_opt)
info['var5'] = list(intercepts_opt)
info['var6'] = list(r_values_opt)
info['var6'] = list(std_errs_opt)
np.savetxt('linreg_full.txt', info, fmt=fmt)

genes_set = []
p_values_set = []
slopes_set = []
intercepts_set = []
r_values_set = []
std_errs_set = []
for gene_id in range(0, len(genes_opt)):

    gene = genes_opt[gene_id]
    p_value = p_values_opt[gene_id]
    slope = slopes_opt[gene_id]
    intercept = intercepts_opt[gene_id]
    r_value = r_values_opt[gene_id]
    std_err = std_errs_opt[gene_id]

    if gene not in genes_set:

        genes_set.append(gene)
        p_values_set.append(p_value)
        slopes_set.append(slope)
        intercepts_set.append(intercept)
        r_values_set.append(r_value)
        std_errs_set.append(std_err)

info = np.zeros(num_top, dtype=[('var1', 'U50'), ('var2', float), ('var3', float), ('var4', float), ('var5', float), ('var6', float)])
fmt = "%s %0.18e %0.18e %0.18e %0.18e %0.18e"
info['var1'] = genes_set[0:num_top]
info['var2'] = p_values_set[0:num_top]
info['var3'] = slopes_set[0:num_top]
info['var4'] = intercepts_set[0:num_top]
info['var5'] = r_values_set[0:num_top]
info['var6'] = std_errs_set[0:num_top]
np.savetxt('linreg_top.txt', info, fmt=fmt)

genes_match = []
for gene in genes_set[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('linreg_match.txt', genes_match, fmt="%s")

print('top: ' + str(len(genes_match)))