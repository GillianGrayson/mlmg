import numpy as np
import scipy.stats as stats
from geo import *
from Infrastructure.load import *
from Infrastructure.save import *
from linreg_mult.routines import *
from method import *
from Infrastructure.file_system import *
from geo import *

method = Method.anova
gd_type = GeneDataType.mean

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

num_top = 100
shift_atr = 5

attributes = get_attributes(config)

min_atr = min(attributes)
max_atr = max(attributes)

min_atr = int(min_atr / shift_atr) * shift_atr
max_atr = (int(max_atr / shift_atr) + 1) * shift_atr
age_dict = {}
for age_id in range(0, len(attributes)):
    age = attributes[age_id]
    key = int((age - min_atr) / shift_atr)
    if key in age_dict:
        age_dict[key].append(age_id)
    else:
        age_dict[key] = [age_id]

gene_names, gene_vals = get_gene_data(config, gd_type)

pvals = []
for id in range(0, len(gene_names)):

    vals = gene_vals[id]

    vals_dict = {}
    for key_age in age_dict:
        vals_dict[key_age] = list(np.asarray(vals)[age_dict[key_age]])

    anova_mean = stats.f_oneway(*vals_dict.values())
    pvals.append(anova_mean.pvalue)

order = np.argsort(pvals)
genes_opt = list(np.array(gene_names)[order])
pvals_opt = list(np.array(pvals)[order])

fn = method.value + '_genes_' + gd_type.value + geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, genes_opt, pvals_opt)

if db_type is DataBaseType.GSE40279:

    table = get_table(config)
    genes_match = []
    for gene in genes_opt[0:num_top]:
        if gene in table:
            genes_match.append(gene)

    fn = method.value + '_match_genes_' + gd_type.value + geo_type.value + '.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)

    print('top: ' + str(len(genes_match)))
