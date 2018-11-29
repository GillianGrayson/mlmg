import statsmodels.api as sm
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from annotations.cpg import *
from annotations.gene import get_dict_cpg_gene
from infrastructure.load.cpg_data import load_dict_cpg_data, get_dict_cpg_gene
from infrastructure.save.features import save_features
from method.clustering.order import *
from infrastructure.load.top import load_top_dict
from method.linreg_mult.routines import *
import pandas as pd


def clock_linreg_mult(config_lvl_2, config_lvl_1):
    attributes = get_attributes(config_lvl_2)
    dict_cpg_data = load_dict_cpg_data(config_lvl_2)
    approved_cpgs = get_approved_cpgs(config_lvl_2)
    dict_cpg_gene = get_dict_cpg_gene(config_lvl_2)

    num_bootstrap_runs = 100

    test_size = math.floor(len(attributes) * config_lvl_2.test_part)
    train_size = len(attributes) - test_size

    if config_lvl_1.method is Method.custom:
        df = pd.read_excel('custom_cpgs.xlsx')
        cpg_names_tmp = list(df['names'])
    else:
        keys = ['cpg'] + get_method_metrics(config_lvl_1.method)
        top_dict = load_top_dict(config_lvl_1, keys)
        cpg_names_tmp = top_dict['cpg']

    cpg_names = []
    cpg_values = []
    for cpg in cpg_names_tmp:
        if cpg in approved_cpgs:
            cpg_names.append(cpg)
            cpg_values.append(dict_cpg_data[cpg])

    keys = ['cpg', 'gene'] + get_method_order_metrics(config_lvl_2.method)
    metrics_dict = {}
    for key in keys:
        metrics_dict[key] = []

    for num_in_top in range(1, min(train_size, len(cpg_names_tmp))):
        print('num_in_top: ' + str(num_in_top))

        metrics_dict['cpg'].append(cpg_names[num_in_top - 1])
        metrics_dict['gene'].append(';'.join(dict_cpg_gene[cpg_names[num_in_top - 1]]))
        metrics_dict['count'].append(num_in_top)

        curr_values = cpg_values[0:num_in_top]

        reg_res = linreg_mult(attributes, curr_values)
        metrics_dict['R2'].append(reg_res.rsquared)

        rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
        indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()

        r_test = 0.0
        evs_test = 0.0
        mae_test = 0.0

        bootstrap_id = 0
        for train_index, test_index in rs.split(indexes):

            X_train = np.array(curr_values).T[train_index].T.tolist()
            X_test = np.array(curr_values).T[test_index].tolist()
            y_train = list(np.array(attributes)[train_index])
            y_test = list(np.array(attributes)[test_index])

            model = linreg_mult(y_train, X_train)

            y_test_pred = model.get_prediction(X_test).predicted_mean
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

        metrics_dict['r_test'].append(r_test)
        metrics_dict['evs_test'].append(evs_test)
        metrics_dict['mae_test'].append(mae_test)

    df = pd.DataFrame(metrics_dict)
    fn = get_result_path(config_lvl_2, 'clock_method(' + config_lvl_1.method.value + ').xlsx')
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

