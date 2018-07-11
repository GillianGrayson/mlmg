from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *

num_top = 100

method = Method.spearman
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

rhos = []
for id in range(0, len(cpgs)):
    curr_vals = vals[id]

    rho, pval = stats.spearmanr(attributes, curr_vals)

    rhos.append(rho)

order = np.argsort(list(map(abs, rhos)))[::-1]
cpgs_sorted = list(np.array(cpgs)[order])
rhos_sorted = list(np.array(rhos)[order])

genes_sorted = []
rhos_genes = []
for id in range(0, len(cpgs_sorted)):
    cpg = cpgs_sorted[id]
    rho = rhos_sorted[id]
    genes = dict_cpg_gene.get(cpg)
    for gene in genes:
        genes_sorted.append(gene)
        rhos_genes.append(rho)

cpgs_sorted = cpgs_sorted[0:num_top]
rhos_sorted = rhos_sorted[0:num_top]

genes_sorted = genes_sorted[0:num_top]
rhos_genes = rhos_genes[0:num_top]

fn = 'cpg/' + method.value + '/' + method.value + '_cpgs.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, cpgs_sorted, rhos_sorted)

fn = 'cpg/' + method.value + '/' + method.value + '_genes_cpg.txt'
fn = get_result_path(fs_type, db_type, fn)
save_params(fn, genes_sorted, rhos_genes)

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