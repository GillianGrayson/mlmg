import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dicts

type = FSType.local_big
print_rate = 1000
num_opt = 1000
num_top = 100

fn = 'table.txt'
full_path = get_path(type, fn)
file = open(full_path)
table = file.read().splitlines()

with open(full_path) as f:
    for line in f:
        table.append(str(line))

suffix = '_ws'

dict_cpg_gene = get_dicts(type)

fn = 'ages.txt'
ages = []
full_path = get_path(type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

fn = 'GSE40279_average_beta.txt'
full_path = get_path(type, fn)

num_lines = 0
genes_rho_dict = {}
genes_abs_rho_dict = {}
genes_pval_dict = {}

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

            for gene in genes:
                if gene in genes_rho_dict:
                    genes_rho_dict[gene] += rho
                else:
                    genes_rho_dict[gene] = rho

                if gene in genes_abs_rho_dict:
                    genes_abs_rho_dict[gene] += abs(rho)
                else:
                    genes_abs_rho_dict[gene] = abs(rho)

                if gene in genes_pval_dict:
                    genes_pval_dict[gene] += pval
                else:
                    genes_pval_dict[gene] = pval

            num_lines += 1

            if num_lines % print_rate == 0:
                print('num_lines: ' + str(num_lines))

rho_list_keys = list(genes_rho_dict.keys())
rho_list_vals = list(genes_rho_dict.values())
rho_order = np.argsort(list(map(abs, rho_list_vals)))[::-1]
rho_list_keys_s = list(np.array(rho_list_keys)[rho_order])
rho_list_vals_s = list(np.array(rho_list_vals)[rho_order])
info = np.zeros(len(rho_list_keys_s), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(rho_list_keys_s)
info['var2'] = list(rho_list_vals_s)
np.savetxt('rho.txt', info, fmt=fmt)
rho_top = 0
for gene in rho_list_keys_s[0:num_top]:
    if gene in table:
        rho_top += 1
print('rho_top: ' + str(rho_top))

abs_rho_list_keys = list(genes_abs_rho_dict.keys())
abs_rho_list_vals = list(genes_abs_rho_dict.values())
abs_rho_order = np.argsort(abs_rho_list_vals)[::-1]
abs_rho_list_keys_s = list(np.array(abs_rho_list_keys)[abs_rho_order])
abs_rho_list_vals_s = list(np.array(abs_rho_list_vals)[abs_rho_order])
info = np.zeros(len(abs_rho_list_keys_s), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(abs_rho_list_keys_s)
info['var2'] = list(abs_rho_list_vals_s)
np.savetxt('abs_rho.txt', info, fmt=fmt)
abs_rho_top = 0
for gene in abs_rho_list_keys_s[0:num_top]:
    if gene in table:
        abs_rho_top += 1
print('abs_rho_top: ' + str(abs_rho_top))

pval_list_keys = list(genes_pval_dict.keys())
pval_list_vals = list(genes_pval_dict.values())
pval_order = np.argsort(pval_list_vals)
pval_list_keys_s = list(np.array(pval_list_keys)[pval_order])
pval_list_vals_s = list(np.array(pval_list_vals)[pval_order])
info = np.zeros(len(pval_list_keys_s), dtype=[('var1', 'U50'), ('var2', float)])
fmt = "%s %0.18e"
info['var1'] = list(pval_list_keys_s)
info['var2'] = list(pval_list_vals_s)
np.savetxt('pval_rho.txt', info, fmt=fmt)
pval_top = 0
for gene in pval_list_keys_s[0:num_top]:
    if gene in table:
        pval_top += 1
print('pval_top: ' + str(pval_top))

order = np.argsort(list(map(abs, rhos)))[::-1]
rhos_opt = list(np.array(rhos)[order[0:num_opt]])
cpgs_opt = list(np.array(cpgs)[order[0:num_opt]])
pvals_opt = list(np.array(pvals)[order[0:num_opt]])
genes_opt = []
for cpg in cpgs_opt:
    genes = dict_cpg_gene.get(cpg)
    sum = ''
    for gene in genes:
        sum += str(gene)
    genes_opt.append(sum)


info = np.zeros(len(cpgs_opt), dtype=[('var1', 'U50'), ('var2', 'U50'), ('var3', float), ('var4', float)])
fmt = "%s %s %0.18e %0.18e"
info['var1'] = list(cpgs_opt)
info['var2'] = list(genes_opt)
info['var3'] = list(rhos_opt)
info['var4'] = list(pvals_opt)
np.savetxt('spearman_full.txt', info, fmt=fmt)

genes_set = []
rhos_set = []
for gene_id in range(0, len(genes_opt)):
    gene = genes_opt[gene_id]
    rho = rhos_opt[gene_id]
    if gene not in genes_set:
        genes_set.append(gene)
        rhos_set.append(rho)
np.savetxt('spearman_genes.txt', genes_set[0:num_top], fmt="%s")
np.savetxt('spearman_rhos.txt', rhos_set[0:num_top], fmt="%0.18e")

genes_match = []
for gene in genes_set[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('spearman_match.txt', genes_match, fmt="%s")

print('top: ' + str(len(genes_match)))