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
from sklearn.model_selection import ShuffleSplit
from enet.routines import *
from Infrastructure.load import *
from Infrastructure.save import *
from linreg_mult.routines import *

num_top_cpgs = 100
train_size = 482
test_size = 174
num_bootstrap_runs = 500

method = Method.enet
fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

attributes = get_attributes(config)

cpgs_top, vals_top = get_top_cpg_data(config, method, num_top_cpgs)

counts, R2s = R2_from_count(vals_top, attributes)

fn = method.value + '_R2s_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_R2(fn, counts, R2s)

metrics_names, metrics_vals = validation_metrics(vals_top, attributes, test_size, train_size, num_bootstrap_runs)

fn = method.value +'_metrics_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, metrics_names, metrics_vals)

print(linreg_mult_with_const(attributes, vals_top).summary())

