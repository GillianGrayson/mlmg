from enum import Enum


class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    linreg_with_rejection = 'linreg_with_rejection'
    linreg_bend = 'linreg_bend'
    linreg_dispersion = 'linreg_dispersion'
    linreg_variance = 'linreg_variance'
    anova = 'anova'
    spearman = 'spearman'
    manova = 'manova'
    random_forest = 'random_forest'
    k_means = 'k_means'
    mean_shift = 'mean_shift'
    linreg_mult = 'linreg_mult'
    match = 'match'

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

def get_method_metrics(method):
    metrics = []
    if method is Method.linreg:
        metrics = ['cluster_mean_shift',
                   'cluster_affinity_prop',
                   'r_value',
                   'p_values',
                   'slope',
                   'intercept']
    elif method is Method.linreg_with_rejection:
        metrics = ['cluster_mean_shift',
                   'cluster_affinity_prop',
                   'r_value',
                   'p_values',
                   'slope',
                   'intercept']
    elif method is Method.linreg_with_rejection:
        metrics = ['cluster_mean_shift',
                   'cluster_affinity_prop',
                   'r_value',
                   'p_values',
                   'slope',
                   'intercept']
    elif method is Method.linreg_bend:
        metrics = ['angles',
                   'slope_left',
                   'intercept_left',
                   'r_value_left',
                   'p_value_left',
                   'std_err_left',
                   'slope_right',
                   'intercept_right',
                   'r_value_right',
                   'p_value_right',
                   'std_err_right']
    elif method is Method.linreg_dispersion:
        metrics = ['slope_left',
                   'intercept_left',
                   'r_value_left',
                   'p_value_left',
                   'std_err_left',
                   'slope_right',
                   'intercept_right',
                   'r_value_right',
                   'p_value_right',
                   'std_err_right']
    elif method is Method.linreg_variance:
        metrics = ['cluster_mean_shift',
                   'cluster_affinity_prop',
                   'r_value',
                   'p_values',
                   'slope',
                   'intercept',
                   'r_value_diff',
                   'p_values_diff',
                   'slope_diff',
                   'intercept_diff']
    elif method is Method.manova:
        metrics = ['pval']

    return metrics
