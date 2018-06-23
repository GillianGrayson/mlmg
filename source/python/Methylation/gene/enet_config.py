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

type = 'mean'

num_folds = 10

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type)

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

regr = ElasticNetCV(cv=num_folds)
elastic_net_X = np.array(vals_passed).T.tolist()
regr.fit(elastic_net_X, attributes)
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






