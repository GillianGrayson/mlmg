import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
import operator
from dicts import *
import statsmodels.api as sm
from sklearn.linear_model import ElasticNetCV, ElasticNet
from config import *
from sklearn.model_selection import ShuffleSplit
from enet.routines import *
from Infrastructure.load import *
from Infrastructure.save import *

train_size = 482
test_size = 174
num_top = 100
num_bootstrap_runs = 100

host_name = socket.gethostname()
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

dict_cpg_gene, dict_cpg_map = get_dicts(config)

fn = 'enet_params_cpg.txt'
path = get_param_path(fs_type, db_type, fn)
params = np.loadtxt(path, dtype='U50')
alpha = float(params[0][1])
l1_ratio = float(params[1][1])

attributes = get_attributes(config)

cpgs_passed, vals_passed = get_cpg_data(config)

rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()
enet_X = np.array(vals_passed).T.tolist()

bootstrap_id = 0
cpg_top_dict = {}
for train_index, test_index in rs.split(indexes):
    print('bootstrap_id: ' + str(bootstrap_id))

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)

    enet_X_train = list(np.array(enet_X)[train_index])
    enet_X_test = list(np.array(enet_X)[test_index])
    enet_y_train = list(np.array(attributes)[train_index])
    enet_y_test = list(np.array(attributes)[test_index])

    enet = enet.fit(enet_X_train, enet_y_train)
    coef = enet.coef_

    order = np.argsort(list(map(abs, coef)))[::-1]
    coef_sorted = list(np.array(coef)[order])
    cpg_sorted = list(np.array(cpgs_passed)[order])
    coef_top = coef_sorted[0:num_top]
    cpg_top = cpg_sorted[0:num_top]
    gene_sorted = []
    for id in range(0, len(cpg_sorted)):
        cpg = cpg_sorted[id]
        genes = dict_cpg_gene.get(cpg)
        for gene in genes:
            if gene not in gene_sorted:
                gene_sorted.append(gene)
    gene_top = gene_sorted[0:num_top]

    for top_id in range(0, num_top):
        cpg = cpg_top[top_id]
        if cpg in cpg_top_dict:
            cpg_top_dict[cpg] += 1
        else:
            cpg_top_dict[cpg] = 1

    bootstrap_id += 1

cpgs = list(cpg_top_dict.keys())
counts = list(cpg_top_dict.values())
order = np.argsort(list(map(abs, counts)))[::-1]
cpgs_sorted = list(np.array(cpgs)[order])
counts_sorted = list(np.array(counts)[order])
genes_sorted = []
counts_genes = []
for id in range(0, len(cpgs_sorted)):
    cpg = cpgs_sorted[id]
    count = counts_sorted[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        genes_sorted.append(gene)
        counts_genes.append(count)

fn = 'enet_cpgs.txt'
fn = get_result_path(fs_type, db_type, fn)
save_enet_top(fn, cpgs_sorted, counts_sorted)

fn = 'enet_genes_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_enet_top(fn, genes_sorted, counts_genes)

if db_type is DataBaseType.GSE40279:
    table = get_table(config)
    genes_match = []
    for gene in list(set(genes_sorted))[0:num_top]:
        if gene in table:
            genes_match.append(gene)
    fn = 'enet_match_genes_cpg.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)
    print('top: ' + str(len(genes_match)))

