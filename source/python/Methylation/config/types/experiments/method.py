from enum import Enum
import math


class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    linreg_with_rejection = 'linreg_with_rejection'
    linreg_bend = 'linreg_bend'
    linreg_dispersion = 'linreg_dispersion'
    linreg_variance = 'linreg_variance'
    linreg_variance_ols = 'linreg_variance_ols'
    linreg_ols = 'linreg_ols'
    linreg_ols_wo_outliers = 'linreg_ols_wo_outliers'
    anova = 'anova'
    anova_statsmodels = 'anova_statsmodels'
    spearman = 'spearman'
    manova = 'manova'
    random_forest = 'random_forest'
    k_means = 'k_means'
    mean_shift = 'mean_shift'
    linreg_mult = 'linreg_mult'
    match = 'match'
    gender_specific = 'gender_specific'
    moment = 'moment'
    classification = 'classification'
    custom = 'custom'

def get_top_fn(method, params_dict):
    fn = 'top.txt'
    if method is Method.linreg_bend:
        fn = 'top'
        for key in params_dict:
            fn += '_' + key + '(' + str(format(params_dict[key], '0.4e')) + ')'
        fn += '.txt'
    elif method is Method.linreg_dispersion:
        fn = 'top'
        for key in params_dict:
            fn += '_' + key + '(' + str(format(params_dict[key], '0.4e')) + ')'
        fn += '.txt'
    return fn

def get_method_metrics(method, is_clustering=False):
    metrics = []
    clustering_metrics = [
        'cluster_mean_shift',
        'cluster_affinity_prop'
    ]
    if method is Method.linreg:
        metrics = [
            'r_value',
            'p_values',
            'slope',
            'intercept'
        ]
        if is_clustering:
            metrics = metrics + clustering_metrics
    elif method is Method.linreg_ols:
        metrics = [
            'R2',
            'intercept',
            'slope',
            'intercepts_std',
            'slopes_std',
            'intercepts_p_values',
            'slopes_p_values'
        ]
        if is_clustering:
            metrics += clustering_metrics
    elif method is Method.linreg_variance_ols:
        metrics = [
            'R2',
            'intercept',
            'slope',
            'intercepts_std',
            'slopes_std',
            'intercepts_p_values',
            'slopes_p_values',
            'R2_var',
            'intercept_var',
            'slope_var',
            'intercepts_std_var',
            'slopes_std_var',
            'intercepts_p_values_var',
            'slopes_p_values_var'
        ]
        if is_clustering:
            metrics += clustering_metrics
    elif method is Method.linreg_with_rejection:
        metrics = [
            'r_value',
            'p_values',
            'slope',
            'intercept'
        ]
        if is_clustering:
            metrics = metrics + clustering_metrics
    elif method is Method.linreg_with_rejection:
        metrics = [
            'r_value',
            'p_values',
            'slope',
            'intercept'
        ]
        if is_clustering:
            metrics = metrics + clustering_metrics
    elif method is Method.linreg_bend:
        metrics = [
            'angles',
            'slope_left',
            'intercept_left',
            'r_value_left',
            'p_value_left',
            'std_err_left',
            'slope_right',
            'intercept_right',
            'r_value_right',
            'p_value_right',
            'std_err_right'
        ]
    elif method is Method.linreg_dispersion:
        metrics = [
            'slope_left',
            'intercept_left',
            'r_value_left',
            'p_value_left',
            'std_err_left',
            'slope_right',
            'intercept_right',
            'r_value_right',
            'p_value_right',
            'std_err_right'
        ]
    elif method is Method.linreg_variance:
        metrics = [
            'r_value',
            'p_values',
            'slope',
            'intercept',
            'r_value_diff',
            'p_values_diff',
            'slope_diff',
            'intercept_diff'
        ]
        if is_clustering:
            metrics = metrics + clustering_metrics
    elif method is Method.manova:
        metrics = ['pval']

    return metrics

def get_method_order_metrics(method):
    metrics = []
    if method is Method.linreg_ols:
        metrics = [
            'names',
            'areas',
            'areas_normed',
            'variance',
            'slope_intersection'
        ]
    elif method is Method.linreg_ols_wo_outliers:
        metrics = [
            'names',
            'areas',
            'areas_normed',
            'variance',
            'slope_intersection'
        ]
    elif method is Method.classification:
        metrics = [
            'names',
            'metric'
        ]
    elif method is Method.manova:
        metrics = [
            'names',
            'gender',
            'age',
            'gender_x_age'
        ]
    elif method is Method.linreg_mult:
        metrics = [
            'count',
            'R2',
            'r_test',
            'evs_test',
            'mae_test'
        ]

    return metrics

def get_method_main_metric(method):
    metric = ''
    if method is Method.linreg:
        metric = 'r_value'
    elif method is Method.manova:
        metric = 'pval'
    return metric

def metric_processing(method, metric):
    if method is Method.linreg:
        return float(metric)
    elif method is Method.manova:
        return -math.log10(float(metric))
    else:
        return float(metric)

def get_method_suffix(params_dict):
    fn = ''
    params_keys = list(params_dict.keys())
    if len(params_keys) > 0:
        params_keys.sort()
        for key in params_keys:
            fn += '_' + key + '(' + str(params_dict[key]) + ')'
    return fn
