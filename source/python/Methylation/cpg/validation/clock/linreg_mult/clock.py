import statsmodels.api as sm
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from annotations.cpg import *
from annotations.gene import get_dict_cpg_gene
from infrastructure.load.cpg_data import load_dict_cpg_data, get_dict_cpg_gene
from infrastructure.save.features import save_features
from method.clustering.order import *
from infrastructure.load.top import load_top_dict
from method.linreg_mult.routines import *
from cpg.validation.clock.clock import *
from itertools import combinations
from config.types.auxiliary import *
import pandas as pd


def build_clock(clock):
    endog_data = clock.endog_data
    endog_name = clock.endog_names
    exog_data = clock.exog_data
    exog_names = clock.exog_names
    metrics_dict = clock.metrics_dict
    train_size = clock.train_size
    test_size = clock.test_size
    num_bootstrap_runs = clock.num_bootstrap_runs
    exog_num = clock.exog_num
    exog_num_comb = clock.exog_num_comb

    endog_dict = {endog_name: endog_data}
    endog_df = pd.DataFrame(endog_dict)

    if exog_num_comb > exog_num:
        exog_num_comb = exog_num

    exog_ids_all = combinations(list(range(0, exog_num)), exog_num_comb)

    R2_best = 0
    r_best = 0
    evs_best = 0
    mae_best = max(endog_data)

    num_comb = 0
    for exog_ids in exog_ids_all:
        num_comb += 1

        exog_dict = {}
        for exog_id in list(exog_ids):
            exog_dict[exog_names[exog_id]] = exog_data[exog_id]
        exog_df = pd.DataFrame(exog_dict)
        exog_df['const'] = 1
        exog_arg_list = ['const'] + exog_names

        reg_res = sm.OLS(endog=endog_df[endog_name], exog=exog_df[exog_arg_list]).fit()

        metrics_dict['summary'].append(reg_res.summary())

        rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
        indexes = np.linspace(0, len(endog_data) - 1, len(endog_data), dtype=int).tolist()

        R2 = reg_res.rsquared
        r_test = 0.0
        evs_test = 0.0
        mae_test = 0.0

        bootstrap_id = 0
        for train_index, test_index in rs.split(indexes):

            endog_train_dict = {endog_name: list(np.array(endog_data)[train_index])}
            endog_train_df = pd.DataFrame(endog_train_dict)

            exog_train_dict = {}
            for exog_id in list(exog_ids):
                exog_train_dict[exog_names[exog_id]] = np.array(exog_data[exog_id]).T[train_index].T.tolist()
            exog_train_df = pd.DataFrame(exog_train_dict)
            exog_train_df['const'] = 1

            y_test = list(np.array(endog_data)[test_index])

            exog_test_dict = {}
            for exog_id in list(exog_ids):
                exog_test_dict[exog_names[exog_id]] = np.array(exog_data[exog_id]).T[test_index].T.tolist()
            exog_test_df = pd.DataFrame(exog_test_dict)
            exog_test_df['const'] = 1

            model = sm.OLS(endog=endog_train_df[endog_name], exog=exog_train_df[exog_arg_list]).fit()

            y_test_pred = model.get_prediction(exog=exog_test_df[exog_arg_list]).predicted_mean
            slope, intercept, r_value, p_value, std_err = stats.linregress(y_test_pred, y_test)
            r_test += r_value
            evs = metrics.explained_variance_score(y_test, list(y_test_pred))
            mae = metrics.mean_absolute_error(y_test, list(y_test_pred))
            evs_test += evs
            mae_test += mae

            bootstrap_id += 1

        r_test /= float(num_bootstrap_runs)
        evs_test /= float(num_bootstrap_runs)
        mae_test /= float(num_bootstrap_runs)

        if mae_test < mae_best:
            R2_best = R2
            r_best = r_test
            evs_best = evs_test
            mae_best = mae_test

        print('num_comb: ' + str(num_comb))

    metrics_dict['R2'].append(R2_best)
    metrics_dict['r_test'].append(r_best)
    metrics_dict['evs_test'].append(evs_best)
    metrics_dict['mae_test'].append(mae_best)

def clock_linreg_mult(config_lvl_2, config_lvl_1):
    attributes = get_attributes(config_lvl_2)
    dict_cpg_data = load_dict_cpg_data(config_lvl_2)
    approved_cpgs = get_approved_cpgs(config_lvl_2)
    dict_cpg_gene = get_dict_cpg_gene(config_lvl_2)

    test_size = math.floor(len(attributes) * config_lvl_2.test_part)
    train_size = len(attributes) - test_size

    if config_lvl_1.method is Method.custom:
        df = pd.read_excel('custom_cpgs.xlsx')
        cpg_names_tmp = list(df['names'])
    else:
        keys = ['cpg'] + get_method_metrics(config_lvl_1.method)

        suffix = get_method_suffix(config_lvl_1.method_params)
        fn = 'top' + suffix + '.txt'

        top_dict = load_top_dict(config_lvl_1, keys, fn=fn)
        cpg_names_tmp = top_dict['cpg']

    cpg_names = []
    cpg_values = []
    for cpg in cpg_names_tmp:
        if cpg in approved_cpgs:
            cpg_names.append(cpg)
            cpg_values.append(dict_cpg_data[cpg])
            if len(cpg_names) >= train_size:
                break

    keys = ['cpg', 'gene'] + get_method_order_metrics(config_lvl_2.method) + ['summary']
    metrics_dict = {}
    for key in keys:
        metrics_dict[key] = []

    if config_lvl_2.method_params is not None:
        exog_type = config_lvl_2.method_params['exog_type']
        exog_num =  min(train_size, len(cpg_names), config_lvl_2.method_params['exog_num'])
        exog_num_comb = min(train_size, len(cpg_names), config_lvl_2.method_params['exog_num_comb'])
    else:
        exog_num = min(train_size, len(cpg_names))
        exog_type = ClockExogType.all
        exog_num_comb = exog_num

    suffix = 'type(' + exog_type.value \
             + ')_num(' + str(exog_num) \
             + ')_num_comb(' + str(exog_num_comb) + ')'

    if exog_type is ClockExogType.all:

        for exog_id in range(0, exog_num):

            print('exog_id: ' + str(exog_id))

            metrics_dict['cpg'].append(cpg_names[exog_id])
            metrics_dict['gene'].append(';'.join(dict_cpg_gene[cpg_names[exog_id]]))
            metrics_dict['count'].append(exog_id + 1)

            clock = Clock(endog_data=attributes,
                          endog_names='age',
                          exog_data=cpg_values[0:exog_id + 1],
                          exog_names=cpg_names[0:exog_id + 1],
                          metrics_dict=metrics_dict,
                          train_size=train_size,
                          test_size=test_size,
                          exog_num=exog_id + 1,
                          exog_num_comb=exog_num_comb)

            build_clock(clock)

    elif exog_type is ClockExogType.single:
        print('exog_num: ' + str(exog_num))
        print('exog_num_comb: ' + str(exog_num_comb))

        metrics_dict['cpg'].append(exog_num_comb)
        metrics_dict['gene'].append(exog_num_comb)
        metrics_dict['count'].append(exog_num_comb)

        clock = Clock(endog_data=attributes,
                      endog_names='age',
                      exog_data=cpg_values[0:exog_num],
                      exog_names=cpg_names[0:exog_num],
                      metrics_dict=metrics_dict,
                      train_size=train_size,
                      test_size=test_size,
                      exog_num=exog_num,
                      exog_num_comb=exog_num_comb)

        build_clock(clock)

    elif exog_type is ClockExogType.slide:
        print('exog_num: ' + str(exog_num))
        print('exog_num_comb: ' + str(exog_num_comb))

        for exog_id in range(0, exog_num, exog_num_comb):

            print('exog_id: ' + str(exog_id))

            metrics_dict['cpg'].append(exog_id)
            metrics_dict['gene'].append(exog_id)
            metrics_dict['count'].append(exog_num_comb)

            clock = Clock(endog_data=attributes,
                          endog_names='age',
                          exog_data=cpg_values[exog_id : exog_id + exog_num_comb],
                          exog_names=cpg_names[exog_id : exog_id + exog_num_comb],
                          metrics_dict=metrics_dict,
                          train_size=train_size,
                          test_size=test_size,
                          exog_num=exog_num_comb,
                          exog_num_comb=exog_num_comb)

            build_clock(clock)

    df = pd.DataFrame(metrics_dict)
    fn = get_result_path(config_lvl_2, 'clock_method(' + config_lvl_1.method.value + ')_' + suffix + '.xlsx')
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

