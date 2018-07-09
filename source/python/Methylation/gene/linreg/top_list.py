import numpy as np
import scipy.stats as stats
from Infrastructure.file_system import *
from geo import *
from Infrastructure.load import *
from Infrastructure.save import *
from linreg_mult.routines import *
from method import *

method = Method.linreg
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

sort_type = 1 # 0 - pval, 1 - rval
num_top = 100

attributes = get_attributes(config)

genes, vals = get_gene_data(config, gd_type)

p_values = []
r_values = []
slopes = []
intercepts = []
for id in range(0, len(genes)):

    val = vals[id]

    slope, intercept, r_value, p_value, std_err = stats.linregress(val, attributes)
    r_values.append(r_value)
    p_values.append(p_value)
    slopes.append(slope)
    intercepts.append(intercept)

order_mean= []
if sort_type == 0:
    order_mean = np.argsort(p_values)
else:
    order_mean = np.argsort(list(map(abs, r_values)))[::-1]
p_values_opt = list(np.array(p_values)[order_mean])
r_values_opt = list(np.array(r_values)[order_mean])
slopes_opt = list(np.array(slopes)[order_mean])
intercepts_opt = list(np.array(intercepts)[order_mean])
genes_opt = list(np.array(genes)[order_mean])

fn = 'gene/' + method.value + '/' +method.value + '_genes_' + gd_type.value + geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_linreg_top(fn, genes_opt, p_values_opt, r_values_opt, slopes_opt, intercepts_opt)

if db_type is DataBaseType.GSE40279:

    table = get_table(config)
    genes_match = []
    for gene in genes_opt[0:num_top]:
        if gene in table:
            genes_match.append(gene)

    fn = 'gene/' + method.value + '/' +method.value + '_match_genes_' + gd_type.value + geo_type.value + '.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)

    print('top: ' + str(len(genes_match)))
