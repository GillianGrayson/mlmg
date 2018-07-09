import math
import numpy as np
import scipy.stats as stats
from Infrastructure.file_system import *
from geo import *
from Infrastructure.load import *
from Infrastructure.save import *
from linreg_mult.routines import *
from method import *
from Infrastructure.file_system import *
from geo import *

method = Method.spearman
gd_type = GeneDataType.mean_der_normed

host_name = socket.gethostname()
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.islands_shores
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

num_top = 100

attributes = get_attributes(config)

gene_names, gene_vals = get_gene_data(config, gd_type)

gene_rhos = []

for id in range(0, len(gene_names)):

    vals = gene_vals[id]
    rhos, pvals = stats.spearmanr(attributes, vals)
    if math.isnan(rhos):
        rhos = 0.0
    gene_rhos.append(rhos)

order = np.argsort(list(map(abs, gene_rhos)))[::-1]
rhos_opt = list(np.array(gene_rhos)[order])
genes_opt = list(np.array(gene_names)[order])

fn = 'gene/' + method.value + '/' +method.value + '_genes_' + gd_type.value + geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, genes_opt, rhos_opt)

if db_type is DataBaseType.GSE40279:

    table = get_table(config)
    genes_match = []
    for gene in genes_opt[0:num_top]:
        if gene in table:
            genes_match.append(gene)

    fn = 'gene/' + method.value + '/' + method.value + '_match_genes_' + gd_type.value + geo_type.value + '.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)

    print('top: ' + str(len(genes_match)))