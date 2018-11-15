from annotations.gene import get_dict_cpg_gene
from infrastructure.path.path import get_path
from config.config import *
from infrastructure.load.routines import line_proc
import numpy as np
import os.path
import pickle


def get_non_inc_cpgs(config):
    cpg_non_inc = []
    if config.data_base is DataBase.GSE40279:
        cpg_non_inc = []
    elif config.data_base is DataBase.GSE52588:
        pval_lim = 0.05
        pval_part = 0.75

        cpgs, pvals = load_cpg_pval_data(config)
        cpg_non_inc = []
        for id in range(0, len(cpgs)):
            curr_cpg = cpgs[id]
            curr_pvals = pvals[id]

            num_big_pvals = 0
            for pval in curr_pvals:
                if pval > pval_lim:
                    num_big_pvals += 1

            if float(num_big_pvals) / float(len(curr_pvals)) > pval_part:
                cpg_non_inc.append(curr_cpg)
    return cpg_non_inc

def load_dict_cpg_data(config):
    indexes = config.indexes

    fn_txt = get_path(config, 'average_beta.txt')
    fn_pkl = get_path(config, 'dict_cpg_data.pkl')

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        dict_cpg_data = pickle.load(f)
        f.close()
    else:
        dict_cpg_data = {}

        f = open(fn_txt)
        for skip_id in range(0, config.num_skip_lines):
             f.readline()

        cpg_non_inc = get_non_inc_cpgs(config)

        num_lines = 0
        for line in f:

            col_vals = line_proc(config, line)

            is_none = False
            if config.miss_tag in col_vals:
                is_none = True

            if not is_none:
                cpg = col_vals[0]
                vals = list(map(float, col_vals[1::]))

                if cpg not in cpg_non_inc:
                    dict_cpg_data[cpg] = vals

            num_lines += 1
            if num_lines % config.print_rate == 0:
                print('num_lines: ' + str(num_lines))
        f.close()

        f = open(fn_pkl, 'wb')
        pickle.dump(dict_cpg_data, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    for cpg, data in dict_cpg_data.items():
        dict_cpg_data[cpg] = list(np.array(data)[indexes])

    return dict_cpg_data

def load_cpg_data_raw(config):
    indexes = config.indexes

    fn = 'average_beta.txt'
    full_path = get_path(config, fn)
    f = open(full_path)
    for skip_id in range(0, config.num_skip_lines):
        skip_line = f.readline()

    num_lines = 0
    cpgs_passed = []
    vals_passed = []

    for line in f:
        col_vals = line_proc(config, line)
        cpg = col_vals[0]
        vals = col_vals[1::]
        vals[-1] = vals[-1].rstrip()

        vals = list(np.array(vals)[indexes])

        vals_passed.append(vals)
        cpgs_passed.append(cpg)

        num_lines += 1
        if num_lines % config.print_rate == 0:
            print('num_lines: ' + str(num_lines))

    f.close()

    return cpgs_passed, vals_passed

def load_cpg_pval_data(config):
    indexes = config.indexes

    fn = 'raw_data.txt'
    fn = get_path(config, fn)
    f = open(fn)
    for skip_id in range(0, config.num_skip_lines):
        skip_line = f.readline()

    num_lines = 0
    cpgs_passed = []
    vals_passed = []
    for line in f:
        col_vals = line.split('\t')
        cpg = col_vals[0].replace('"', '')
        vals = list(map(float, col_vals[1::]))
        pvals = vals[2::3]
        pvals = list(np.array(pvals)[indexes])

        cpgs_passed.append(cpg)
        vals_passed.append(pvals)

        num_lines += 1
        if num_lines % config.print_rate == 0:
            print('num_lines: ' + str(num_lines))

    f.close()

    return cpgs_passed, vals_passed

