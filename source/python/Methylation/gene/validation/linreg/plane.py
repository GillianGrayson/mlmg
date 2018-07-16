import scipy.stats as stats

from config import *
from infrastructure.load import *
from infrastructure.save import *
from method.method import *

num_top_genes = 100

method = Method.enet
val_method = Validation.linreg
gd_type_order = GeneDataType.mean
gd_type_vals = GeneDataType.mean_der_normed
host_name = socket.gethostname()
geo_type = GeoType.islands_shores
fs_type = FSType.local_big
if host_name == 'MSI':
    fs_type = FSType.local_msi
elif host_name == 'DESKTOP-K9VO2TI':
    fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

attributes = get_attributes(config)

gene_names = get_top_gene_names(config, gd_type_order, method, num_top_genes)
gene_vals_main= get_top_gene_vals(config, gd_type_order, gene_names)
gene_vals_aux = get_top_gene_vals(config, gd_type_vals, gene_names)

p_values_main = []
r_values_main = []
p_values_aux = []
r_values_aux = []
for id in range(0, len(gene_names)):
    vals_main = gene_vals_main[id]
    slope, intercept, r_value, p_value, std_err = stats.linregress(vals_main, attributes)
    p_values_main.append(p_value)
    r_values_main.append(r_value)

    vals_aux = gene_vals_aux[id]
    slope, intercept, r_value, p_value, std_err = stats.linregress(vals_aux, attributes)
    p_values_aux.append(p_value)
    r_values_aux.append(r_value)

fn = 'gene/' + val_method.value + '/' + method.value + '_plane(' + gd_type_order.value + ')' +  geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_features(fn, gene_names, [r_values_main, r_values_aux] )

