from sklearn.linear_model import ElasticNetCV
from config import *
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