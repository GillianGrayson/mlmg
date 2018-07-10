from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *

num_top = 100
shift_atr = 5

method = Method.anova
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

cpgs, vals = get_cpg_data(config)

pvals = []
for id in range(0, len(cpgs)):
    curr_vals = vals[id]

    curr_beta_dict = {}
    for key_age in age_dict:
        curr_beta_dict[key_age] = list(np.asarray(curr_vals)[age_dict[key_age]])

    anova_res = stats.f_oneway(*curr_beta_dict.values())
    pvals.append(anova_res.pvalue)

order = np.argsort(pvals)
cpgs_sorted = list(np.array(cpgs)[order])
pvals_sorted = list(np.array(pvals)[order])
genes_sorted = []
pvals_genes = []
for id in range(0, len(cpgs_sorted)):
    cpg = cpgs_sorted[id]
    pval = pvals_sorted[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        genes_sorted.append(gene)
        pvals_genes.append(pval)

cpgs_sorted = cpgs_sorted[0:num_top]
pvals_sorted = pvals_sorted[0:num_top]

genes_sorted = genes_sorted[0:num_top]
pvals_genes = pvals_genes[0:num_top]

fn = 'cpg/' + method.value + '/' + method.value + '_cpgs.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, cpgs_sorted, pvals_sorted)

fn = 'cpg/' + method.value + '/' + method.value + '_genes_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, genes_sorted, pvals_genes)

if db_type is DataBaseType.GSE40279:
    table = get_table(config)
    genes_match = []
    for gene in list(set(genes_sorted))[0:num_top]:
        if gene in table:
            genes_match.append(gene)
    fn = 'cpg/' + method.value + '/' + method.value + '_match_genes_cpg.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)
    print('top: ' + str(len(genes_match)))