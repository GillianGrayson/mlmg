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

def mult_linreg_with_const(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

def mult_linreg(y, x):
    results = sm.OLS(y, np.array(x).T).fit()
    return results

top = 100
bootstrap = 100

num_top_cpgs = 100

train_size = 482
test_size = 174

num_bootstrap_runs = 100

print_rate = 10000

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

cpgs_top = []
fn = 'enet_cpgs_top(' + str(top) + ')_bootstrap(' + str(bootstrap) + ').txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for line in f:
    cpg = line.split(' ')[0]
    cpgs_top.append(cpg)

cpgs_top = cpgs_top[0:num_top_cpgs]

fn = db_type.value + '_average_beta.txt'
full_path = get_full_path(fs_type, db_type, fn)
f = open(full_path)
for skip_id in range(0, config.num_skip_lines):
    skip_line = f.readline()

num_lines = 0
dict_top = {}

for line in f:

    col_vals = config.line_proc(line)
    cpg = col_vals[0]
    vals = list(map(float, col_vals[1::]))

    if cpg in cpgs_top:
        dict_top[cpg] = vals

    num_lines += 1
    if num_lines % print_rate == 0:
        print('num_lines: ' + str(num_lines))

vals_top = []
for cpg in cpgs_top:
    vals = dict_top.get(cpg)
    vals_top.append(vals)

R2s = []
nums = []

for num_opt in range(0, len(vals_top)):
    reg_res = mult_linreg(attributes, vals_top[0:num_opt + 1])
    nums.append(num_opt + 1)
    R2s.append(reg_res.rsquared)

info = np.zeros(len(nums), dtype=[('var1', int), ('var2', float)])
fmt = "%d %0.18e"
info['var1'] = nums
info['var2'] = R2s
np.savetxt('R2s.txt', info, fmt=fmt)

rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()

bootstrap_id = 0
r_avg_test = 0.0
std_err_avg_test = 0.0
r_avg_train = 0.0
std_err_avg_train = 0.0
cpg_top_dict = {}
for train_index, test_index in rs.split(indexes):
    print('bootstrap_id: ' + str(bootstrap_id))

    X_train = np.array(vals_top).T[train_index].T.tolist()
    X_test = np.array(vals_top).T[test_index].tolist()
    y_train = list(np.array(attributes)[train_index])
    y_test = list(np.array(attributes)[test_index])

    model = mult_linreg(y_train, X_train)

    y_test_pred = model.get_prediction(X_test).predicted_mean
    slope, intercept, r_value, p_value, std_err = stats.linregress(y_test_pred, y_test)
    r_avg_test += r_value
    std_err_avg_test += std_err

    y_train_pred = model.get_prediction(list(np.array(X_train).T)).predicted_mean
    slope, intercept, r_value, p_value, std_err = stats.linregress(y_train_pred, y_train)
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


