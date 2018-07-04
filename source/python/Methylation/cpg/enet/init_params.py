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
from enet.routines import *
from Infrastructure.load import *
from Infrastructure.save import *

num_folds = 10

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

attributes = get_attributes(fs_type, db_type)

cpgs, vals = get_cpg_data(config)

param_names, param_values = get_enet_params(attributes, vals, num_folds)

fn = 'enet_params_cpg.txt'
fn = get_param_path(fs_type, db_type, fn)
save_params(fn, param_names, param_values)

