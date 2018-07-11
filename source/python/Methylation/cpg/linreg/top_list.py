from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *

num_top = 100

method = Method.linreg
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

dict_cpg_gene = get_dict_cpg_gene(config)

attributes = get_attributes(config)

cpgs, vals = get_cpg_data(config)

slopes = []
intercepts = []
rvals = []
pvals = []
for id in range(0, len(cpgs)):
    curr_vals = vals[id]

    slope, intercept, r_value, p_value, std_err = stats.linregress(curr_vals, attributes)

    slopes.append(slope)
    intercepts.append(intercept)
    rvals.append(r_value)
    pvals.append(p_value)

order = np.argsort(pvals)
cpgs_sorted = list(np.array(cpgs)[order])
pvals_sorted = list(np.array(pvals)[order])
slopes_sorted = list(np.array(slopes)[order])
intercepts_sorted = list(np.array(intercepts)[order])
rvals_sorted = list(np.array(rvals)[order])

genes_sorted = []
pvals_genes = []
slopes_genes = []
intercepts_genes = []
rvals_genes = []
for id in range(0, len(cpgs_sorted)):
    cpg = cpgs_sorted[id]
    pval = pvals_sorted[id]
    slope = slopes_sorted[id]
    intercept = intercepts_sorted[id]
    rval = rvals_sorted[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        genes_sorted.append(gene)
        pvals_genes.append(pval)
        slopes_genes.append(slope)
        intercepts_genes.append(intercept)
        rvals_genes.append(rval)

cpgs_sorted = cpgs_sorted[0:num_top]
pvals_sorted = pvals_sorted[0:num_top]
slopes_sorted = slopes_sorted[0:num_top]
intercepts_sorted = intercepts_sorted[0:num_top]
rvals_sorted = rvals_sorted[0:num_top]

genes_sorted = genes_sorted[0:num_top]
pvals_genes = pvals_genes[0:num_top]
slopes_genes = slopes_genes[0:num_top]
intercepts_genes = intercepts_genes[0:num_top]
rvals_genes = rvals_genes[0:num_top]

fn = 'cpg/' + method.value + '/' + method.value + '_cpgs.txt'
fn = get_result_path(fs_type, db_type, fn)
save_linreg_top(fn, cpgs_sorted, pvals_sorted, rvals_sorted, slopes_sorted, intercepts_sorted)

fn = 'cpg/' + method.value + '/' + method.value + '_genes_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_linreg_top(fn, genes_sorted, pvals_genes, rvals_genes, slopes_genes, intercepts_genes)

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
