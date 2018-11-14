from infrastructure.load.routines import line_proc
from infrastructure.path.path import *
from config.types.experiments.method import *
import numpy as np
import os.path
import pandas as pd


def load_top_gene_names_by_article(config, fn):
    full_path = get_origin_path(config, fn)
    file = open(full_path)
    table = file.read().splitlines()
    file.close()
    return table

def load_top_gene_data(config, num_top):
    gene_names = load_top_gene_names(config, num_top)
    gene_vals = load_top_gene_vals(config, gene_names)
    return gene_names, gene_vals

def load_top_gene_names(config, num_top):
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    f = open(fn)
    gene_names = []
    for line in f:
        gene = line.split(' ')[0].rstrip()
        gene_names.append(gene)
    gene_names = gene_names[0:num_top]
    return gene_names

def load_top_data(config, num_top, index):
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    f = open(fn)
    gene_names = []
    for line in f:
        gene = line.split(' ')[index].rstrip()
        gene_names.append(gene)
    gene_names = gene_names[0:num_top]
    return gene_names

def load_top_dict(config, keys, fn='top.txt', num_top=1000000):
    if fn.endswith('.txt'):
        params_last_run_dict = load_top_params_last_run_dict(config)
        if params_last_run_dict:
            fn = get_top_fn(config.approach_method, params_last_run_dict)
        fn = get_result_path(config, fn)
        f = open(fn)
        top_dict = {}
        num_lines = 0
        for line in f:
            line = line.rstrip()
            cols = line.split(' ')
            for key_id in range(0, len(keys)):
                key = keys[key_id]
                if key not in top_dict:
                    top_dict[key] = [cols[key_id]]
                else:
                    top_dict[key].append(cols[key_id])
            num_lines += 1
            if num_lines >= num_top:
                print('Top contains more than ' + str(num_top))
                break
        return top_dict
    elif fn.endswith('.xlsx'):
        fn = get_result_path(config, fn)
        df = pd.read_excel(fn)
        top_dict = {}
        for key in keys:
            top_dict[key] = list(df[key])
        return top_dict

def load_top_params_last_run_dict(config):
    fn = 'params_last_run.txt'
    fn = get_result_path(config, fn)
    params_last_run_dict = {}
    if os.path.isfile(fn):
        f = open(fn)
        for line in f:
            cols = line.split(' ')
            key = cols[0]
            val = float(cols[1])
            params_last_run_dict[key] = val
    return params_last_run_dict

def load_top_gene_vals(config, genes_top):
    indexes = config.indexes
    dict_top = {}
    fn = 'gene_data.txt'
    fn = get_data_path(config, fn)
    f = open(fn)
    for line in f:
        col_vals = line.split(' ')
        gene = col_vals[0]
        vals = list(map(float, col_vals[1::]))
        vals = list(np.array(vals)[indexes])
        if gene in genes_top:
            dict_top[gene] = vals
    gene_vals = []
    for gene in genes_top:
        vals = dict_top.get(gene)
        gene_vals.append(vals)
    return gene_vals

def load_top_gene_names_by_cpg(config, method, num_top):
    fn = 'genes_from_cpg.txt'
    fn = get_result_path(config, fn)
    f = open(fn)
    gene_names = []
    for line in f:
        gene = line.split(' ')[0].rstrip()
        gene_names.append(gene)
    gene_names = gene_names[0:num_top]
    return gene_names

def load_top_cpg_data(config, method, num_top):
    indexes = config.indexes
    db_type = config.db_type
    print_rate = config.print_rate
    cpgs_top = []
    fn = 'top.txt'
    fn = get_result_path(config, fn)
    f = open(fn)
    for line in f:
        cpg = line.split(' ')[0].rstrip()
        cpgs_top.append(cpg)

    cpgs_top = cpgs_top[0:num_top]

    fn = db_type.value + '_average_beta.txt'
    path = get_path(config, fn)
    f = open(path)
    for skip_id in range(0, config.num_skip_lines):
        skip_line = f.readline()

    num_lines = 0
    dict_top = {}

    for line in f:

        col_vals = line_proc(config, line)
        cpg = col_vals[0]
        vals = list(map(float, col_vals[1::]))
        vals = list(np.array(vals)[indexes])

        if cpg in cpgs_top:
            dict_top[cpg] = vals

        num_lines += 1
        if num_lines % print_rate == 0:
            print('num_lines: ' + str(num_lines))

    vals_top = []
    for cpg in cpgs_top:
        vals = dict_top.get(cpg)
        vals_top.append(vals)

    return cpgs_top, vals_top

