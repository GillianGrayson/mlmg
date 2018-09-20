from config.config import *
from annotations.bop import *
from infrastructure.load.cpg_data import load_cpg_data
from statsmodels.multivariate.manova import MANOVA
from statsmodels.stats.multitest import multipletests
from infrastructure.save.features import save_features
from infrastructure.load.bop_data import load_bop_cpg_dict
import pandas as pd


def save_top_manova(config, attributes_types, attribute_target, num_top=500, window=3, test=MANOVATest.pillai_bartlett):
    dict_bop_cpgs = load_bop_cpg_dict(config)
    dict_bop_genes = get_dict_bop_genes(config, dict_bop_cpgs)
    cpgs, betas = load_cpg_data(config)

    atr_table = []
    atr_cols = []
    for atr_type in attributes_types:
        if isinstance(atr_type, Attribute):
            atr_table.append(get_attributes(config, atr_type))
        elif isinstance(atr_type, CellPop):
            atr_table.append(get_cell_pop(config, [atr_type]))
        atr_cols.append(atr_type.value)

    num_bops = 0
    bops_passed = []
    bops_pvals = []
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

                table = list(map(list, zip(*table)))
                x = pd.DataFrame(table, columns=cols)
                manova = MANOVA.from_formula(formula, x)
                mv_test_res = manova.mv_test()
                pvals = mv_test_res.results[attribute_target.value]['stat'].values[0:4, 4]
                target_pval = pvals[0]
                if test is MANOVATest.wilks:
                    target_pval = pvals[0]
                elif test is MANOVATest.pillai_bartlett:
                    target_pval = pvals[1]
                elif test is MANOVATest.lawley_hotelling:
                    target_pval = pvals[2]
                elif test is MANOVATest.roy:
                    target_pval = pvals[3]
                pvals_on_bop.append(target_pval)
            min_pval = np.min(pvals_on_bop)
            bops_passed.append(bop)
            bops_pvals.append(min_pval)
        num_bops += 1
        if num_bops % config.print_rate == 0:
            print('num_bops: ' + str(num_bops))

    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(bops_pvals, 0.05, method='fdr_bh')
    order = np.argsort(pvals_corrected)
    bops_opt = list(np.array(bops_passed)[order])[0:num_top]
    pvals_opt = list(np.array(pvals_corrected)[order])[0:num_top]
    genes_opt = []
    genes_from_bop = []
    for bop in bops_opt:
        curr_genes = dict_bop_genes.get(bop)
        genes_str = curr_genes[0]
        for gene_id in range(1, len(curr_genes)):
            genes_str += ';' + curr_genes[gene_id]
        genes_opt.append(genes_str)
        for gene in curr_genes:
            if gene not in genes_from_bop:
                genes_from_bop.append(gene)

    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [bops_opt, genes_opt, pvals_opt])

    config.approach_gd = GeneDataType.from_bop
    config.dt = DataType.gene
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [genes_from_bop])
    config.dt = DataType.cpg
