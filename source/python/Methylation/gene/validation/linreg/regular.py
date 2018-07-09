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
from sklearn import metrics
from Infrastructure.load import *
from Infrastructure.save import *
from linreg_mult.routines import *
from method import *

num_top_genes = 100

method = Method.enet
val_method = Validation.linreg
gd_type = GeneDataType.mean_der_normed
host_name = socket.gethostname()
geo_type = GeoType.any
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

attributes = get_attributes(config)

gene_names, gene_vals = get_top_gene_data(config, gd_type, method, num_top_genes)

p_values = []
r_values = []
slopes = []
intercepts = []
for id in range(0, len(gene_names)):
    vals = gene_vals[id]
    slope, intercept, r_value, p_value, std_err = stats.linregress(vals, attributes)
    r_values.append(r_value)
    p_values.append(p_value)
    slopes.append(slope)
    intercepts.append(intercept)

order = np.argsort(list(map(abs, r_values)))[::-1]
p_values = list(np.array(p_values)[order])
r_values = list(np.array(r_values)[order])
slopes = list(np.array(slopes)[order])
intercepts = list(np.array(intercepts)[order])
gene_names = list(np.array(gene_names)[order])

fn = 'gene/' + val_method.value + '/' + method.value + '_metrics_' + gd_type.value + geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_linreg_top(fn, gene_names, p_values, r_values, slopes, intercepts)

