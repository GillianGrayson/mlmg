from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from scipy import stats
from attributes.conditions import *
from config.config import *
from copy import deepcopy


def save_bend_linreg(config, limit, pval, num_opt=1000):
    config_less = deepcopy(config)
    age_less(config_less, limit)
    atr_l = get_attributes(config_less)
    dict_cpg_data = load_dict_cpg_data(config_less)
    cpg_names_l = list(dict_cpg_data.keys())
    cpg_vals_l = list(dict_cpg_data.values())

    config_more = deepcopy(config)
    age_more(config_more, limit)
    atr_m = get_attributes(config_more)
    dict_cpg_data = load_dict_cpg_data(config_more)
    cpg_names_m = list(dict_cpg_data.keys())
    cpg_vals_m = list(dict_cpg_data.values())

    cpg_gene_dict = get_dict_cpg_gene(config)

    cpgs_passed = []
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

    num_cpgs = 0

    for cpg_id_l in range(0, len(cpg_names_l)):
        cpg_id_m = cpg_names_m.index(cpg_names_l[cpg_id_l])
        vals_l = cpg_vals_l[cpg_id_l]
        vals_m = cpg_vals_m[cpg_id_m]

        slope_l, intercept_l, r_value_l, p_value_l, std_err_l = stats.linregress(atr_l, vals_l)
        slope_m, intercept_m, r_value_m, p_value_m, std_err_m = stats.linregress(atr_m, vals_m)
        angle = abs(slope_l - slope_m)

        if (max(p_value_l, p_value_m) < pval):
            cpgs_passed.append(cpg_names_l[cpg_id_l])

            genes = cpg_gene_dict.get(cpg_names_l[cpg_id_l])
            if len(genes) > 0:
                if genes[0] == '':
                    genes_passed.append('nan')
                else:
                    genes_passed.append(";".join(genes))
            else:
                genes_passed.append('nan')

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

        num_cpgs += 1
        if num_cpgs % config.print_rate == 0:
            print('num_cpgs: ' + str(num_cpgs))

    order = np.argsort(angles)[::-1][0:num_opt]

    cpgs_opt = list(np.array(cpgs_passed)[order])

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

    fn = get_result_path(config, 'bend_' + str(limit) + '.txt')
    save_features(fn, [cpgs_opt,

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
                       std_err_ms_opt])

    raw_config = Config(db=config.db,
                        dt=config.dt,
                        approach=config.approach,
                        scenario=config.scenario,
                        approach_method=config.approach_method,
                        gender=Gender.any)

    cpg_str_list = []
    dict_cpg_data = load_dict_cpg_data(raw_config)
    cpg_name_raw = list(dict_cpg_data.keys())
    cpg_vals_raw = list(dict_cpg_data.values())
    for cpg in cpgs_opt:
        cpg_vals = cpg_vals_raw[cpg_name_raw.index(cpg)]
        curr_cpg_str = cpg
        for id in range(0, len(cpg_vals)):
            curr_cpg_str += (' ' + str(format(cpg_vals[id], '0.8e')))
        cpg_str_list.append(curr_cpg_str)

    fn = get_result_path(config, 'bend_data_' + str(limit) + '.txt')
    np.savetxt(fn, cpg_str_list, fmt="%s")


db = DataBase.GSE40279
dt = DataType.cpg
approach = Approach.bend
scenario = Scenario.approach
approach_method = Method.linreg
gender = Gender.F

limit = 55
pval = 0.00001

config = Config(
    db=db,
    dt=dt,
    approach=approach,
    scenario=scenario,
    approach_method=approach_method,
    gender=gender
)

save_bend_linreg(config, limit, pval)
