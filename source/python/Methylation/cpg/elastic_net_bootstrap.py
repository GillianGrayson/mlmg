import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
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

alpha = 1.564782107693006921e-02
l1_ratio = 5.000000000000000000e-01

print_rate = 10000
num_top = 100
num_bootstrap_runs = 10

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

rs = ShuffleSplit(num_bootstrap_runs, 174, 482)
indexes = np.linspace(0, len(ages) - 1, len(ages), dtype=int).tolist()
enet_X = np.array(vals_passed).T.tolist()

bootstrap_id = 0
r_avg = 0.0
std_err_avg = 0.0
cpg_top_dict = {}
for train_index, test_index in rs.split(indexes):
    print('bootstrap_id: ' + str(bootstrap_id))

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)

    enet_X_train = list(np.array(enet_X)[train_index])
    enet_X_test = list(np.array(enet_X)[test_index])
    enet_y_train = list(np.array(ages)[train_index])
    enet_y_test = list(np.array(ages)[test_index])

    enet = enet.fit(enet_X_train, enet_y_train)
    coef = enet.coef_

    order = np.argsort(list(map(abs, coef)))[::-1]
    cpg_sorted = list(np.array(cpgs_passed)[order])
    cpg_top = cpg_sorted[0:num_top]
    coef_sorted = list(np.array(coef)[order])
    coef_top = coef_sorted[0:num_top]
    gene_sorted = []
    for id in range(0, len(cpg_sorted)):
        cpg = cpg_sorted[id]
        genes = dict_cpg_gene.get(cpg)
        for gene in genes:
            if gene not in gene_sorted:
                gene_sorted.append(gene)
    gene_top = gene_sorted[0:num_top]

    for top_id in range(0, num_top):
        cpg = cpg_top[top_id]
        if cpg in cpg_top_dict:
            cpg_top_dict[cpg] += 1
        else:
            cpg_top_dict[cpg] = 0

    enet_y_pred = enet.predict(enet_X_test).tolist()

    slope, intercept, r_value, p_value, std_err = stats.linregress(enet_y_pred, enet_y_test)
    r_avg += r_value
    std_err_avg += std_err

    bootstrap_id += 1

r_avg /= float(num_bootstrap_runs)
std_err_avg /= float( num_bootstrap_runs)
print('r_avg: ', r_avg)
print('std_err_avg: ', std_err_avg)

sorted_cpg = sorted(cpg_top_dict, key=cpg_top_dict.get, reverse=True)



ololo = 1



