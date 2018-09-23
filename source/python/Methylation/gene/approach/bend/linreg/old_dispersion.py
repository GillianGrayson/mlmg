from method.enet.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.path import get_result_path
from infrastructure.save.features import save_features
from scipy import stats
from attributes.conditions import *
from config.config import *
from copy import deepcopy


def save_old_dispersion(config, limit, young_pval):
    config_less = deepcopy(config)
    age_less(config_less, limit)
    atr_l = get_attributes(config_less)
    g_names_l, g_vals_l = load_gene_data(config_less)

    config_more = deepcopy(config)
    age_more(config_more, limit)
    atr_m = get_attributes(config_more)
    g_names_m, g_vals_m = load_gene_data(config_more)

    genes_passed = []

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
        vals_l = g_vals_l[g_id_l]
        slope_l, intercept_l, r_value_l, p_value_l, std_err_l = stats.linregress(atr_l, vals_l)
        if (p_value_l < young_pval):
            genes_passed.append(g_names_l[g_id_l])

            g_id_m = g_names_m.index(g_names_l[g_id_l])
            vals_m = g_vals_m[g_id_m]
            slope_m, intercept_m, r_value_m, p_value_m, std_err_m = stats.linregress(atr_m, vals_m)

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

    order = np.argsort(p_value_ms)[::-1]

    genes_opt = list(np.array(genes_passed)[order])

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

    fn = get_result_path(config, 'old_dispersion_age(' + str(limit) +')' + '.txt')
    save_features(fn, [genes_opt,

                       slope_ls_opt,
                       intercept_ls_opt,
                       r_value_ls_opt,
                       p_value_ls_opt,
                       std_err_ls_opt,

                       slope_ms_opt,
                       intercept_ms_opt,
                       r_value_ms_opt,
                       p_value_ms_opt,
                       std_err_ms_opt])


db = DataBaseType.GSE40279
dt = DataType.gene
approach = Approach.bend
scenario = Scenario.approach
approach_method = Method.linreg
approach_gd = GeneDataType.mean
genders = [Gender.F, Gender.M, Gender.any]
geo = GeoType.islands_shores

limit = 65
pval = 1e-7

for gender in genders:
    print('gender: ' + str(gender.value))

    config = Config(
        db=db,
        dt=dt,
        approach=approach,
        scenario=scenario,
        approach_method=approach_method,
        approach_gd=approach_gd,
        geo=geo,
        gender=gender
    )

    save_old_dispersion(config, limit, pval)
