import numpy as np
from config import *
from infrastructure.file_system import *
from annotation.regular import *

db_type = DataBaseType.GSE52588
geo_type = GeoType.any
host_name = socket.gethostname()
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
config = Config(fs_type, db_type, geo_type=geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type=geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type=geo_type)

gene_cpg_dict = get_dict_gene_cpg(config)
cpgs, vals = get_cpg_data(config)

fn = 'graph_genes_2'

f = open(fn + '.txt')
target_genes = f.read().splitlines()

str_list = []
for gene in target_genes:
    gene_cpgs = gene_cpg_dict[gene]
    for gene_cpg in gene_cpgs:
        if gene_cpg in cpgs:
            index = cpgs.index(gene_cpg)
            curr_vals = vals[cpgs.index(gene_cpg)]
            curr_str = gene + ' ' + gene_cpg
            for id in range(0, len(curr_vals)):
                curr_str += (' ' + str(format(curr_vals[id], '0.8e')))
            str_list.append(curr_str)

np.savetxt(fn + '_cpgs.txt', str_list, fmt="%s")
