import numpy as np
from config.config import *
from attributes.categorical import get_attributes_dict
from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.path.path import get_result_path
from infrastructure.save.features import save_features
from annotations.gene import get_dict_cpg_gene
from annotations.cpg import *
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import pandas as pd


def save_top_anova_statsmodels(config):
    attributes_dict = get_attributes_dict(config)
    dict_cpg_gene = get_dict_cpg_gene(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpg_names = list(dict_cpg_data.keys())
    cpg_values = list(dict_cpg_data.values())
    approved_cpgs = get_approved_cpgs(config)

    atr_table = []
    atr_cols = []
    for atr_type in config.attributes_types:
        if isinstance(atr_type, Attribute):
            atr = get_attributes(config, atr_type)
            atr = proceed_attributes(atr, atr_type)
            atr_table.append(atr)
        elif isinstance(atr_type, CellPop):
            atr_table.append(get_cell_pop(config, [atr_type]))
        atr_cols.append(atr_type.value)

    cpg_names_passed = []
    p_values = []
    for id in range(0, len(cpg_names)):
        cpg = cpg_names[id]
        if cpg in approved_cpgs:
            val_table = []
            val_cols = []

            cpg_names_passed.append(cpg)
            curr_vals = cpg_values[id]
            val_table.append(curr_vals)
            val_cols.append('cpg')
            table = atr_table + val_table
            cols = atr_cols + val_cols

            formula = val_cols[0]
            formula += ' ~ ' + atr_cols[0]
            for atr_col_id in range(1, len(atr_cols)):
                atr_col = atr_cols[atr_col_id]
                formula += ' + ' + atr_col
            for target_id in range(0, len(config.attribute_target)):
                if isinstance(config.attribute_target[target_id], tuple):
                    formula += ' + ' + '*'.join([x.value for x in config.attribute_target[target_id]])

            table = list(map(list, zip(*table)))
            x = pd.DataFrame(table, columns=cols)
            anova_model = ols(formula, data = x).fit()
            anova_test_res = sm.stats.anova_lm(anova_model, typ=2)
            pvals = []
            for target_id in range(0, len(config.attribute_target)):
                if isinstance(config.attribute_target[target_id], tuple):
                    atr_name = ':'.join([x.value for x in config.attribute_target[target_id]])
                else:
                    atr_name = config.attribute_target[target_id].value
                pvals.append(anova_test_res.T[atr_name].values[3])
            p_values.append(pvals)
            if id % config.print_rate == 0:
                print('cpg_id: ' + str(id))

    p_values = list(map(list, zip(*p_values)))
    p_values_corrected = []
    for target_id in range(0, len(config.attribute_target)):
        reject, pvals_corr, alphacSidak, alphacBonf = multipletests(p_values[target_id], 0.05, method='fdr_bh')
        p_values_corrected.append(pvals_corr)
    order = np.argsort(p_values_corrected[-1])
    cpgs_sorted = list(np.array(cpg_names_passed)[order])
    pvals_sorted = [np.array(x)[order] for x in p_values_corrected]

    target_str = config.attribute_target[0].value
    is_cells = False
    for target_id in range(1, len(config.attribute_target)):
        if isinstance(config.attribute_target[target_id], tuple):
            tmp = []
            is_mult_cells = False
            for x in config.attribute_target[target_id]:
                if isinstance(x, Attribute):
                    tmp.append(x.value)
                elif isinstance(x, CellPop):
                    is_mult_cells = True
            target_str += '_' + '_x_'.join([x for x in tmp])
            if is_mult_cells:
                target_str += '_x_cells'
        elif isinstance(config.attribute_target[target_id], CellPop):
            if not is_cells:
                is_cells = True
                target_str += '_cells'
        else:
            target_str += '_' + config.attribute_target[target_id].value

    types_str = config.attributes_types[0].value
    is_cells = False
    for type_id in range(1, len(config.attributes_types)):
        if isinstance(config.attributes_types[type_id], Attribute):
            types_str += '_' + config.attributes_types[type_id].value
        elif isinstance(config.attributes_types[type_id], CellPop):
            if not is_cells:
                is_cells = True
                types_str += '_cells'

    suffix = '_target(' + target_str + ')_exog(' + types_str + ').txt'

    features = [cpgs_sorted]
    [features.append(list(x)) for x in pvals_sorted]
    fn = 'top' + suffix
    fn = get_result_path(config, fn)
    save_features(fn, features)

    genes_sorted = []
    pvals_genes = []
    for target_id in range(0, len(config.attribute_target)):
        pvals_genes.append([])
    for id in range(0, len(cpgs_sorted)):
        cpg = cpgs_sorted[id]
        p_value = []
        for target_id in range(0, len(config.attribute_target)):
            p_value.append(pvals_sorted[target_id][id])
        if cpg in dict_cpg_gene:
            genes = dict_cpg_gene.get(cpg)
            for gene in genes:
                if gene not in genes_sorted:
                    genes_sorted.append(gene)
                    for target_id in range(0, len(config.attribute_target)):
                        pvals_genes[target_id].append(p_value[target_id])

    features = [genes_sorted]
    [features.append(list(x)) for x in pvals_genes]

    config.gene_data_type = GeneDataType.from_cpg
    config.data_type = DataType.gene
    fn = 'top' + suffix
    fn = get_result_path(config, fn)
    save_features(fn, features)
    config.data_type = DataType.cpg
