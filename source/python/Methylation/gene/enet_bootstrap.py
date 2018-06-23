import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
import operator
from dicts import *
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV, ElasticNet
from config import *
from sklearn.model_selection import ShuffleSplit

def reg_m(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

type = 'mean'

alpha = 8.454092531050155530e-04
l1_ratio = 5.000000000000000000e-01

train_size = 482
test_size = 174

print_rate = 10000
num_top = 100
num_bootstrap_runs = 100

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
attributes = []
full_path = get_full_path(fs_type, db_type, fn)
with open(full_path) as f:
    for line in f:
        attributes.append(int(line))

genes_passed = []
vals_passed = []
fn = 'gene_' + type + geo_type.value + '.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    col_vals = line.split(' ')
    gene = col_vals[0]
    vals = list(map(float, col_vals[1::]))
    genes_passed.append(gene)
    vals_passed.append(vals)

rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()
enet_X = np.array(vals_passed).T.tolist()

bootstrap_id = 0
r_avg_test = 0.0
std_err_avg_test = 0.0
r_avg_train = 0.0
std_err_avg_train = 0.0
gene_top_dict = {}
for train_index, test_index in rs.split(indexes):
    print('bootstrap_id: ' + str(bootstrap_id))

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)

    enet_X_train = list(np.array(enet_X)[train_index])
    enet_X_test = list(np.array(enet_X)[test_index])
    enet_y_train = list(np.array(attributes)[train_index])
    enet_y_test = list(np.array(attributes)[test_index])

    enet = enet.fit(enet_X_train, enet_y_train)
    coef = enet.coef_

    order = np.argsort(list(map(abs, coef)))[::-1]
    coef_sorted = list(np.array(coef)[order])
    gene_sorted = list(np.array(genes_passed)[order])
    coef_top = coef_sorted[0:num_top]
    gene_top = gene_sorted[0:num_top]

    for top_id in range(0, num_top):
        gene = gene_top[top_id]
        if gene in gene_top_dict:
            gene_top_dict[gene] += 1
        else:
            gene_top_dict[gene] = 1

    enet_y_test_pred = enet.predict(enet_X_test).tolist()
    slope, intercept, r_value, p_value, std_err = stats.linregress(enet_y_test_pred, enet_y_test)
    r_avg_test += r_value
    std_err_avg_test += std_err

    enet_y_train_pred = enet.predict(enet_X_train).tolist()
    slope, intercept, r_value, p_value, std_err = stats.linregress(enet_y_train_pred, enet_y_train)
    r_avg_train += r_value
    std_err_avg_train += std_err

    bootstrap_id += 1

r_avg_test /= float(num_bootstrap_runs)
std_err_avg_test /= float(num_bootstrap_runs)
r_avg_train /= float(num_bootstrap_runs)
std_err_avg_train /= float(num_bootstrap_runs)
print('r_avg_test: ', r_avg_test)
print('std_err_avg_test: ', std_err_avg_test)
print('r_avg_train: ', r_avg_train)
print('std_err_avg_train: ', std_err_avg_train)

genes = list(gene_top_dict.keys())
counts = list(gene_top_dict.values())
order = np.argsort(list(map(abs, counts)))[::-1]
genes_sorted = list(np.array(genes)[order])
counts_sorted = list(np.array(counts)[order])

info = np.zeros(len(genes_sorted), dtype=[('var1', 'U50'),  ('var2', int)])
fmt = "%s %d"
info['var1'] = list(genes_sorted)
info['var2'] = list(counts_sorted)
np.savetxt('enet_bootstrap_genes.txt', info, fmt=fmt)

