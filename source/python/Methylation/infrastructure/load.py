import numpy as np
from annotation.regular import *
from method.method import *
from infrastructure.file_system import *

def get_attributes(config):
    db_type = config.db_type
    fs_type = config.fs_type
    fn = 'attribute.txt'
    attributes = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            if db_type is DataBaseType.GSE40279:
                attributes.append(int(line))
            elif db_type is DataBaseType.GSE52588:
                attributes.append(int(line))
    return attributes

def get_cpg_data(config):
    db_type = config.db_type
    fs_type = config.fs_type
    print_rate = config.print_rate

    dict_cpg_gene = get_dict_cpg_gene(config)

    fn = db_type.value + '_average_beta.txt'
    full_path = get_path(fs_type, db_type, fn)
    f = open(full_path)
    for skip_id in range(0, config.num_skip_lines):
        skip_line = f.readline()

    num_lines = 0
    cpgs_passed = []
    vals_passed = []

    for line in f:

        col_vals = config.line_proc(line)

        is_none = False
        if config.miss_tag in col_vals:
            is_none = True

        if not is_none:
            cpg = col_vals[0]
            vals = list(map(float, col_vals[1::]))

            if cpg in dict_cpg_gene:
                vals_passed.append(vals)
                cpgs_passed.append(cpg)

            num_lines += 1
            if num_lines % print_rate == 0:
                print('num_lines: ' + str(num_lines))

    return cpgs_passed, vals_passed

def get_cpg_pval_data(config):
    db_type = config.db_type
    fs_type = config.fs_type

    num_skip_lines_raw = 1

    fn = db_type.value + '_raw_data.txt'
    fn = get_path(fs_type, db_type, fn)
    f = open(fn)
    for skip_id in range(0, num_skip_lines_raw):
        skip_line = f.readline()

    num_lines = 0
    cpgs_passed = []
    vals_passed = []
    for line in f:
        col_vals = line.split('\t')
        cpg = col_vals[0].replace('"', '')
        vals = list(map(float, col_vals[1::]))
        pvals = vals[2::3]

        cpgs_passed.append(cpg)
        vals_passed.append(pvals)

        num_lines += 1
        if num_lines % config.print_rate == 0:
            print('num_lines: ' + str(num_lines))

    return cpgs_passed, vals_passed

def get_gene_data(config, gd_type):
    db_type = config.db_type
    fs_type = config.fs_type
    geo_type = config.geo_type

    genes_passed = []
    vals_passed = []
    fn = 'gene_' + gd_type.value + geo_type.value + '.txt'
    path = get_gene_data_path(fs_type, db_type, gd_type, fn)
    f = open(path)
    for line in f:
        col_vals = line.split(' ')
        gene = col_vals[0]
        vals = list(map(float, col_vals[1::]))
        nans = np.isnan(vals)
        if True not in nans:
            genes_passed.append(gene)
            vals_passed.append(vals)

    return genes_passed, vals_passed

def get_table(config):
    db_type = config.db_type
    fs_type = config.fs_type

    fn = 'table.txt'
    full_path = get_path(fs_type, db_type, fn)
    file = open(full_path)
    table = file.read().splitlines()

    return table

def get_params_dict(config, gd_type, method):
    db_type = config.db_type
    fs_type = config.fs_type
    geo_type = config.geo_type
    fn = method.value + '_params_' + gd_type.value + geo_type.value + '.txt'
    path = get_param_path(fs_type, db_type, fn)
    params = np.loadtxt(path, dtype='U50')

    params_dict = {}
    if method is Method.enet:
        alpha = float(params[0][1])
        l1_ratio = float(params[1][1])
        params_dict['alpha'] = alpha
        params_dict['l1_ratio'] = l1_ratio

    return params_dict

def get_top_gene_data(config, gd_type, method, num_top):
    gene_names = get_top_gene_names(config, gd_type, method, num_top)
    gene_vals = get_top_gene_vals(config, gd_type, gene_names)

    return gene_names, gene_vals

def get_top_gene_names(config, gd_type, method, num_top):
    db_type = config.db_type
    fs_type = config.fs_type
    geo_type = config.geo_type
    fn = ''
    if method is Method.clustering:
        fn = 'gene/' + method.value + '/' + config.clustering_type.value + '_1D_' + gd_type.value + geo_type.value + '.txt'
    else:
        fn = 'gene/' + method.value + '/' + method.value + '_genes_' + gd_type.value + geo_type.value + '.txt'
    fn = get_result_path(fs_type, db_type, fn)
    f = open(fn)
    gene_names = []
    for line in f:
        gene = line.split(' ')[0].rstrip()
        gene_names.append(gene)

        gene_names = gene_names[0:num_top]

    return gene_names

def get_top_cpg_gene_names(config, method, num_top):
    db_type = config.db_type
    fs_type = config.fs_type
    geo_type = config.geo_type
    fn = 'cpg/' + method.value + '/' + method.value + '_genes_cpgs.txt'
    fn = get_result_path(fs_type, db_type, fn)
    f = open(fn)
    gene_names = []
    for line in f:
        gene = line.split(' ')[0].rstrip()
        gene_names.append(gene)

        gene_names = gene_names[0:num_top]

    return gene_names

def get_top_gene_vals(config, gd_type, genes_top):
    db_type = config.db_type
    fs_type = config.fs_type
    geo_type = config.geo_type
    dict_top = {}
    fn = 'gene_' + gd_type.value + geo_type.value + '.txt'
    fn = get_gene_data_path(fs_type, db_type, gd_type, fn)
    f = open(fn)
    for line in f:
        col_vals = line.split(' ')
        gene = col_vals[0]
        vals = list(map(float, col_vals[1::]))
        if gene in genes_top:
            dict_top[gene] = vals

    gene_vals = []
    for gene in genes_top:
        vals = dict_top.get(gene)
        gene_vals.append(vals)

    return gene_vals


def get_top_cpg_data(config, method, num_top):
    db_type = config.db_type
    fs_type = config.fs_type
    print_rate = config.print_rate
    cpgs_top = []
    fn = 'cpg/' + method.value + '/' + method.value + '_cpgs.txt'
    fn = get_result_path(fs_type, db_type, fn)
    f = open(fn)
    for line in f:
        cpg = line.split(' ')[0].rstrip()
        cpgs_top.append(cpg)

    cpgs_top = cpgs_top[0:num_top]

    fn = db_type.value + '_average_beta.txt'
    path = get_path(fs_type, db_type, fn)
    f = open(path)
    for skip_id in range(0, config.num_skip_lines):
        skip_line = f.readline()

    num_lines = 0
    dict_top = {}

    for line in f:

        col_vals = config.line_proc(line)
        cpg = col_vals[0]
        vals = list(map(float, col_vals[1::]))

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

