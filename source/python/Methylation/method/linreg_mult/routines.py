import statsmodels.api as sm
import numpy as np
from sklearn.linear_model import ElasticNetCV, ElasticNet
from config import *
from sklearn.model_selection import ShuffleSplit
from sklearn import metrics
import scipy.stats as stats

def linreg_mult_with_const(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

def linreg_mult(y, x):
    results = sm.OLS(y, np.array(x).T).fit()
    return results
