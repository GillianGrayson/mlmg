from infrastructure.load.gene_data import load_gene_data
from infrastructure.save.features import save_features
from scipy import stats
from attributes.conditions import *
from config.config import *
from copy import deepcopy


def save_top_linreg_bend(config, limit=65, pval=0.001, is_clustering=False):
    config_less = deepcopy(config)
    age_less(config_less, limit)
    atr_l = get_attributes(config_less)
    g_names_l, g_vals_l = load_gene_data(config_less)

    config_more = deepcopy(config)
    age_more(config_more, limit)
    atr_m = get_attributes(config_more)
    g_names_m, g_vals_m = load_gene_data(config_more)

    genes_passed = []

    angles = []

    slope_ls = []
    intercept_ls = []
    r_value_ls = []
    p_value_ls = []
    std_err_ls = []

    slope_ms = []
    intercept_ms = []
    r_value_ms = []
    p_value_ms = []
    std_err_ms = []

    for g_id_l in range(0, len(g_names_l)):
        g_id_m = g_names_m.index(g_names_l[g_id_l])
        vals_l = g_vals_l[g_id_l]
        vals_m = g_vals_m[g_id_m]

        slope_l, intercept_l, r_value_l, p_value_l, std_err_l = stats.linregress(atr_l, vals_l)
        slope_m, intercept_m, r_value_m, p_value_m, std_err_m = stats.linregress(atr_m, vals_m)
        angle = abs(slope_l - slope_m)

        if (max(p_value_l, p_value_m) < pval):
            genes_passed.append(g_names_l[g_id_l])
            angles.append(angle)

            slope_ls.append(slope_l)
            intercept_ls.append(intercept_l)
            r_value_ls.append(r_value_l)
            p_value_ls.append(p_value_l)
            std_err_ls.append(std_err_l)

            slope_ms.append(slope_m)
            intercept_ms.append(intercept_m)
            r_value_ms.append(r_value_m)
            p_value_ms.append(p_value_m)
            std_err_ms.append(std_err_m)

    order = np.argsort(angles)[::-1]
    genes_opt = list(np.array(genes_passed)[order])
    angles_opt = list(np.array(angles)[order])

    slope_ls_opt = list(np.array(slope_ls)[order])
    intercept_ls_opt = list(np.array(intercept_ls)[order])
    r_value_ls_opt = list(np.array(r_value_ls)[order])
    p_value_ls_opt = list(np.array(p_value_ls)[order])
    std_err_ls_opt = list(np.array(std_err_ls)[order])

    slope_ms_opt = list(np.array(slope_ms)[order])
    intercept_ms_opt = list(np.array(intercept_ms)[order])
    r_value_ms_opt = list(np.array(r_value_ms)[order])
    p_value_ms_opt = list(np.array(p_value_ms)[order])
    std_err_ms_opt = list(np.array(std_err_ms)[order])

    features = [
        genes_opt,

        angles_opt,

        slope_ls_opt,
        intercept_ls_opt,
        r_value_ls_opt,
        p_value_ls_opt,
        std_err_ls_opt,

        slope_ms_opt,
        intercept_ms_opt,
        r_value_ms_opt,
        p_value_ms_opt,
        std_err_ms_opt
    ]
    fn = get_result_path(config, 'top_limit(' + str(format(limit, '0.4e')) + ')_pval(' + str(format(pval, '0.4e')) + ')' + '.txt')
    save_features(fn, features)

    fn = get_result_path(config, 'params_last_run.txt')
    params_strs = ['limit' + ' ' + str(format(limit, '0.4e')),
                   'pval' + ' ' + str(format(pval, '0.4e'))]
    np.savetxt(fn, params_strs, fmt="%s")