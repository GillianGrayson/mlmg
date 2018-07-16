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

def R2_from_count(vals, attributes):
    counts = []
    R2s = []
    for num_opt in range(0, len(vals)):
        reg_res = linreg_mult(attributes, vals[0:num_opt + 1])
        counts.append(num_opt + 1)
        R2s.append(reg_res.rsquared)

    return counts, R2s

def validation_metrics(vals, attributes, test_size, train_size, num_bootstrap_runs):

    rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
    indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()

    bootstrap_id = 0

    r_avg_test = 0.0
    std_err_avg_test = 0.0
    r_avg_train = 0.0
    std_err_avg_train = 0.0

    evs_avg_test = 0.0
    mae_avg_test = 0.0
    evs_avg_train = 0.0
    mae_avg_train = 0.0

    for train_index, test_index in rs.split(indexes):
        print('bootstrap_id: ' + str(bootstrap_id))

        X_train = np.array(vals).T[train_index].T.tolist()
        X_test = np.array(vals).T[test_index].tolist()
        y_train = list(np.array(attributes)[train_index])
        y_test = list(np.array(attributes)[test_index])

        model = linreg_mult(y_train, X_train)

        y_test_pred = model.get_prediction(X_test).predicted_mean
        slope, intercept, r_value, p_value, std_err = stats.linregress(y_test_pred, y_test)
        r_avg_test += r_value
        std_err_avg_test += std_err

        evs = metrics.explained_variance_score(y_test, list(y_test_pred))
        mae = metrics.mean_absolute_error(y_test, list(y_test_pred))
        evs_avg_test += evs
        mae_avg_test += mae

        y_train_pred = model.get_prediction(list(np.array(X_train).T)).predicted_mean
        slope, intercept, r_value, p_value, std_err = stats.linregress(y_train_pred, y_train)
        r_avg_train += r_value
        std_err_avg_train += std_err

        evs = metrics.explained_variance_score(y_train, list(y_train_pred))
        mae = metrics.mean_absolute_error(y_train, list(y_train_pred))
        evs_avg_train += evs
        mae_avg_train += mae

        bootstrap_id += 1

    r_avg_test /= float(num_bootstrap_runs)
    std_err_avg_test /= float(num_bootstrap_runs)
    r_avg_train /= float(num_bootstrap_runs)
    std_err_avg_train /= float(num_bootstrap_runs)
    evs_avg_test /= float(num_bootstrap_runs)
    mae_avg_test /= float(num_bootstrap_runs)
    evs_avg_train /= float(num_bootstrap_runs)
    mae_avg_train /= float(num_bootstrap_runs)

    metrics_names = ['r_avg_test', 'std_err_avg_test', 'r_avg_train', 'std_err_avg_train', 'evs_avg_test',
                    'mae_avg_test', 'evs_avg_train', 'mae_avg_train']
    metrics_vals = [r_avg_test, std_err_avg_test, r_avg_train, std_err_avg_train, evs_avg_test, mae_avg_test,
                    evs_avg_train, mae_avg_train]

    return metrics_names, metrics_vals
