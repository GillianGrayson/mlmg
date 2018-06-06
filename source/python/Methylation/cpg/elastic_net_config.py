import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import *
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV
from config import *

print_rate = 10000
num_top = 100
num_folds = 10
num_bootstrap_runs = 500

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type)

dict_cpg_gene, dict_cpg_map = get_dicts(fs_type, db_type, geo_type)

fn = 'attribute.txt'
ages = []
full_path = get_full_path(fs_type, db_type, fn)
with open(full_path) as f:
    for line in f:
        ages.append(int(line))

fn = db_type.value + '_average_beta.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for skip_id in range(0, config.num_skip_lines):
    skip_line = f.readline()

num_lines = 0
cpgs_passed = []
vals_passed = []

for line in f:

    col_vals = line.split('\t')
    cpg = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    if cpg in dict_cpg_gene:
        vals_passed.append(vals)
        cpgs_passed.append(cpg)

    num_lines += 1
    if num_lines % print_rate == 0:
        print('num_lines: ' + str(num_lines))

regr = ElasticNetCV(cv=num_folds)
elastic_net_X = np.array(vals_passed).T.tolist()
regr.fit(elastic_net_X, ages)
coef = regr.coef_
alpha = regr.alpha_
l1_ratio = regr.l1_ratio_
param_names = ['alpha',  'l1_ratio']
param_values = [alpha, l1_ratio]
info = np.zeros(2, dtype=[('var1', 'U50'), ('var2', 'float')])
fmt = "%s %0.18e"
info['var1'] = param_names
info['var2'] = param_values
np.savetxt('elastic_net_params.txt', info, fmt=fmt)

# Local saving: only cpg
order = np.argsort(list(map(abs, coef)))[::-1]
cpg_sorted = list(np.array(cpgs_passed)[order])
cpgs_top = cpg_sorted[0:num_top]
coef_sorted = list(np.array(coef)[order])
coef_top = coef_sorted[0:num_top]
genes_str_top = []

for cpg in cpgs_top:
    genes = dict_cpg_gene.get(cpg)
    genes_str = ''
    for gene in genes:
        genes_str += (str(gene) + ';')
    genes_str = genes_str[:-1]
    genes_str_top.append(genes_str)

info = np.zeros(num_top, dtype=[('var1', 'U50'), ('var2', 'U50'), ('var3', 'float')])
fmt = "%s %s %0.18e"
info['var1'] = cpgs_top
info['var2'] = genes_str_top
info['var3'] = coef_top
np.savetxt('cpg_top.txt', info, fmt=fmt)

# Local saving: only gene
genes_top = []
coef_genes_top = []
for id in range(0, len(cpg_sorted)):
    cpg = cpg_sorted[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        if gene not in genes_top:
            genes_top.append(gene)
            coef_genes_top.append(coef_sorted[id])

info = np.zeros(num_top, dtype=[('var1', 'U50'), ('var2', 'float')])
fmt = "%s %0.18e"
info['var1'] = genes_top[0:num_top]
info['var2'] = coef_genes_top[0:num_top]
np.savetxt('gene_top.txt', info, fmt=fmt)

# Checking genes with table in article 2015
fn = 'table.txt'
full_path = get_full_path(fs_type, db_type, fn)
file = open(full_path)
table = file.read().splitlines()
genes_match = []
for gene in genes_top[0:num_top]:
    if gene in table:
        genes_match.append(gene)
np.savetxt('match.txt', genes_match, fmt="%s")

print('top: ' + str(len(genes_match)))





