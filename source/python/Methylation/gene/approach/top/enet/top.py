from sklearn.linear_model import ElasticNet
from sklearn.model_selection import ShuffleSplit
from method.enet.routines import *
from infrastructure.load.params import load_params_dict
from infrastructure.load.attributes import get_attributes
from infrastructure.load.gene_data import load_gene_data
from infrastructure.file_system import get_result_path
from infrastructure.save.features import save_features


def save_top_enet(config, num_bootstrap_runs=100, num_top=500):

    params_dict = load_params_dict(config)
    alpha = params_dict.get('alpha')
    l1_ratio = params_dict.get('l1_ratio')

    attributes = get_attributes(config)

    genes_passed, vals_passed = load_gene_data(config)

    test_size = int(len(attributes) * config.test_part)
    train_size = len(attributes) - test_size
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

    fn = get_result_path(config, 'top.txt')
    save_features(fn, [genes_sorted, counts_sorted])

