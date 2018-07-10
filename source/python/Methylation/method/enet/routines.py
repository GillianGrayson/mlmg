import pathlib
from infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import *
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV
from config import *

def get_enet_params(x, y, num_folds):
    regr = ElasticNetCV(cv=num_folds)
    elastic_net_X = np.array(x).T.tolist()
    regr.fit(elastic_net_X, y)
    alpha = regr.alpha_
    l1_ratio = regr.l1_ratio_
    names = ['alpha', 'l1_ratio']
    values = [alpha, l1_ratio]
    return names, values