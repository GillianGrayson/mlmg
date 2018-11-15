from config.config import *
from config.types.auxiliary import *
from annotations.bop import *
from config.types.attributes.attribute import Attribute
from config.types.attributes.cell_pop import CellPop
from infrastructure.load.cpg_data import load_dict_cpg_data
from statsmodels.multivariate.manova import MANOVA
from statsmodels.stats.multitest import multipletests
from infrastructure.save.features import save_features
from sklearn.cluster import MeanShift, estimate_bandwidth, AffinityPropagation
from method.clustering.order import *
import pandas as pd


def save_top_manova(config, window=3, test=MANOVATest.pillai_bartlett):
    dict_bop_cpgs = get_dict_bop_cpgs(config)
    dict_bop_genes = get_dict_bop_genes(config)
    dict_cpg_data = load_dict_cpg_data(config)
    cpgs = list(dict_cpg_data.keys())
    betas = list(dict_cpg_data.values())

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

    num_bops = 0
    bops_passed = []
    p_values = []
    for bop in dict_bop_cpgs:
        curr_cpgs = dict_bop_cpgs.get(bop)
        cpgs_passed = []
        for cpg in curr_cpgs:
            if cpg in cpgs:
                cpgs_passed.append(cpg)
        if len(cpgs_passed) > 2:
            pvals_on_bop = []
            for win_id in range(0, len(cpgs_passed) - 2):
                val_table = []
                val_cols = []
                for cpg_id in range(0, window):
                    cpg = cpgs_passed[win_id + cpg_id]
                    beta = betas[cpgs.index(cpg)]
                    val_table.append(beta)
                    val_cols.append('cpg_'+str(cpg_id))
                table = atr_table + val_table
                cols = atr_cols + val_cols

                formula = val_cols[0]
                for val_col_id in range(1, len(val_cols)):
                    val_col = val_cols[val_col_id]
                    formula += ' + ' + val_col
                formula += ' ~ ' + atr_cols[0]
                for atr_col_id in range(1, len(atr_cols)):
                    atr_col = atr_cols[atr_col_id]
                    formula += ' + ' + atr_col
                for target_id in range(0, len(config.attribute_target)):
                    if isinstance(config.attribute_target[target_id],tuple):
                        formula += ' + ' + '*'.join([x.value for x in config.attribute_target[target_id]])

                table = list(map(list, zip(*table)))
                x = pd.DataFrame(table, columns=cols)
                manova = MANOVA.from_formula(formula, x)
                mv_test_res = manova.mv_test()
                pvals = []
                target_pval = []
                for target_id in range(0, len(config.attribute_target)):
                    if isinstance(config.attribute_target[target_id],tuple):
                        atr_name = ':'.join([x.value for x in config.attribute_target[target_id]])
                    else:
                        atr_name = config.attribute_target[target_id].value
                    pvals.append(mv_test_res.results[atr_name]['stat'].values[0:4, 4])
                    if test is MANOVATest.wilks:
                        target_pval.append(pvals[target_id][0])
                    elif test is MANOVATest.pillai_bartlett:
                        target_pval.append(pvals[target_id][1])
                    elif test is MANOVATest.lawley_hotelling:
                        target_pval.append(pvals[target_id][2])
                    elif test is MANOVATest.roy:
                        target_pval.append(pvals[target_id][3])
                pvals_on_bop.append(target_pval)
            pvals_on_bop = list(map(list, zip(*pvals_on_bop)))
            argmin_pval = int(np.argmin(pvals_on_bop[-1]))
            min_pvals = []
            for target_id in range(0, len(config.attribute_target)):
                min_pvals.append(pvals_on_bop[target_id][argmin_pval])
            bops_passed.append(bop)
            p_values.append(min_pvals)
        num_bops += 1
        if num_bops % config.print_rate == 0:
            print('num_bops: ' + str(num_bops))

    p_values = list(map(list, zip(*p_values)))
    p_values_corrected = []
    for target_id in range(0, len(config.attribute_target)):
        reject, pvals_corr, alphacSidak, alphacBonf = multipletests(p_values[target_id], 0.05, method='fdr_bh')
        p_values_corrected.append(pvals_corr)
    order = np.argsort(p_values_corrected[-1])
    bops_sorted = list(np.array(bops_passed)[order])
    p_values_sorted = [np.array(x)[order] for x in p_values_corrected]

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

    clusters_mean_shift = []
    clusters_affinity_prop = []
    features = [bops_sorted]
    [features.append(list(x)) for x in p_values_sorted]
    if config.is_clustering:
        metrics_sorted_np = np.asarray(list(map(np.log10, p_values_sorted))).reshape(-1, 1)
        bandwidth = estimate_bandwidth(metrics_sorted_np)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(metrics_sorted_np)
        labels_mean_shift = list(ms.labels_)
        clusters_mean_shift = clustering_order(labels_mean_shift)
        af = AffinityPropagation().fit(metrics_sorted_np)
        labels_affinity_propagation = list(af.labels_)
        clusters_affinity_prop = clustering_order(labels_affinity_propagation)
        features = features + [
            clusters_mean_shift,
            clusters_affinity_prop
        ]
        fn = 'top_with_clustering' + suffix
    else:
        fn = 'top' + suffix
    fn = get_result_path(config, fn)
    save_features(fn, features)

    genes_sorted = []
    p_values_genes = []
    clusters_mean_shift_genes = []
    clusters_affinity_prop_genes = []
    for id in range(0, len(bops_sorted)):
        bop = bops_sorted[id]
        p_value = p_values_sorted[-1][id]
        if config.is_clustering:
            cluster_mean_shift = clusters_mean_shift[id]
            cluster_affinity_prop = clusters_affinity_prop[id]
        else:
            cluster_mean_shift = []
            cluster_affinity_prop = []
        if bop in dict_bop_genes:
            genes = dict_bop_genes.get(bop)
            for gene in genes:
                if gene not in genes_sorted:
                    genes_sorted.append(gene)
                    p_values_genes.append(p_value)
                    if config.is_clustering:
                        clusters_mean_shift_genes.append(cluster_mean_shift)
                        clusters_affinity_prop_genes.append(cluster_affinity_prop)

    config.data_type = DataType.gene
    gene_data_type = config.gene_data_type
    geo_type = config.geo_type
    config.gene_data_type = GeneDataType.from_bop
    config.geo_type = GeoType.from_bop

    features = [
        genes_sorted,
        p_values_genes,
    ]
    if config.is_clustering:
        fn = 'top_with_clustering' + suffix
        features = features + [
            clusters_mean_shift_genes,
            clusters_affinity_prop_genes
        ]
    else:
        fn = 'top' + suffix

    fn = get_result_path(config, fn)
    save_features(fn, features)

    config.gene_data_type = gene_data_type
    config.geo_type = geo_type
    config.data_type = DataType.bop
