from sklearn.linear_model import ElasticNet
from sklearn.model_selection import ShuffleSplit

from infrastructure.load import *
from infrastructure.save import *
from method.enet.routines import *
from method.method import *

method = Method.enet
gd_type = GeneDataType.mean_der_normed

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
geo_type = GeoType.islands_shores
config = Config(fs_type, db_type, geo_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type, geo_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type, geo_type)

params_dict = get_params_dict(config, gd_type, method)
alpha = params_dict.get('alpha')
l1_ratio = params_dict.get('l1_ratio')

attributes = get_attributes(config)

genes_passed, vals_passed = get_gene_data(config, gd_type)

rs = ShuffleSplit(num_bootstrap_runs, test_size, train_size)
indexes = np.linspace(0, len(attributes) - 1, len(attributes), dtype=int).tolist()
enet_X = np.array(vals_passed).T.tolist()

bootstrap_id = 0
gene_top_dict = {}
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
    gene_sorted = list(np.array(genes_passed)[order])
    coef_top = coef_sorted[0:num_top]
    gene_top = gene_sorted[0:num_top]

    for top_id in range(0, num_top):
        gene = gene_top[top_id]
        if gene in gene_top_dict:
            gene_top_dict[gene] += 1
        else:
            gene_top_dict[gene] = 1

    bootstrap_id += 1

genes = list(gene_top_dict.keys())
counts = list(gene_top_dict.values())
order = np.argsort(list(map(abs, counts)))[::-1]
genes_sorted = list(np.array(genes)[order])
counts_sorted = list(np.array(counts)[order])

fn = 'gene/' + method.value + '/' +method.value + '_genes_' + gd_type.value + geo_type.value + '.txt'
fn = get_result_path(fs_type, db_type, fn)
save_enet_top(fn, genes_sorted, counts_sorted)

if db_type is DataBaseType.GSE40279:
    table = get_table(config)
    genes_match = []
    for gene in genes_sorted[0:num_top]:
        if gene in table:
            genes_match.append(gene)

    fn = 'gene/' + method.value + '/' +method.value + '_match_genes_' + gd_type.value + geo_type.value + '.txt'
    fn = get_result_path(fs_type, db_type, fn)
    save_names(fn, genes_match)

    print('top: ' + str(len(genes_match)))

