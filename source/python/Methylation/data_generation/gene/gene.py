import numpy as np
from config.config import *
from infrastructure.load.gene_data import get_raw_dict
from infrastructure.path import get_gene_data_path
from config.types import *

def save_gene_data(config):

    gene_raw_dict = get_raw_dict(config)

    gene_mean_str_list = []
    gene_std_str_list = []
    gene_mean_der_str_list = []
    gene_mean_der_normed_str_list = []
    num_genes = 0
    for gene in gene_raw_dict:

        mean_list = []
        std_list = []
        mean_der_list = []
        mean_der_normed_list = []
        for curr_list in gene_raw_dict[gene]:
            mean_list.append(np.mean(curr_list))
            std_list.append(np.std(curr_list))

            if len(curr_list) > 1:
                curr_sum = 0.0
                for elem_id in range(0, len(curr_list) - 1):
                    curr_sum += abs(curr_list[elem_id + 1] - curr_list[elem_id])
                curr_sum /= float(len(curr_list) - 1)
                mean_der_list.append(curr_sum)
                tmp_sum = np.sum(curr_list)
                norm_cpg = tmp_sum / float(len(curr_list))
                normed = curr_sum / norm_cpg
                mean_der_normed_list.append(normed)
            else:
                mean_der_list.append(0.0)
                mean_der_normed_list.append(0.0)

        curr_mean_str = gene
        curr_std_str = gene
        curr_mean_der_str = gene
        curr_mean_der_normed_str = gene
        for id in range(0, len(mean_list)):
            curr_mean_str += (' ' + str(format(mean_list[id], '0.8e')))
            curr_std_str += (' ' + str(format(std_list[id], '0.8e')))
            curr_mean_der_str += (' ' + str(format(mean_der_list[id], '0.8e')))
            curr_mean_der_normed_str += (' ' + str(format(mean_der_normed_list[id], '0.8e')))

        gene_mean_str_list.append(curr_mean_str)
        gene_std_str_list.append(curr_std_str)
        gene_mean_der_str_list.append(curr_mean_der_str)
        gene_mean_der_normed_str_list.append(curr_mean_der_normed_str)

        num_genes += 1
        if num_genes % config.print_rate == 0:
            print('num_genes: ' + str(num_genes))

    # mean
    config.approach_gd = GeneDataType.mean
    config.validation_gd = GeneDataType.mean
    fn = get_gene_data_path(config, 'gene_data.txt')
    np.savetxt(fn, gene_mean_str_list, fmt="%s")
    # std
    config.approach_gd = GeneDataType.std
    config.validation_gd = GeneDataType.std
    fn = get_gene_data_path(config, 'gene_data.txt')
    np.savetxt(fn, gene_std_str_list, fmt="%s")
    # mean_der
    config.approach_gd = GeneDataType.mean_der
    config.validation_gd = GeneDataType.mean_der
    fn = get_gene_data_path(config, 'gene_data.txt')
    np.savetxt(fn, gene_mean_der_str_list, fmt="%s")
    # mean_der_normed
    config.approach_gd = GeneDataType.mean_der_normed
    config.validation_gd = GeneDataType.mean_der_normed
    fn = get_gene_data_path(config, 'gene_data.txt')
    np.savetxt(fn, gene_mean_der_normed_str_list, fmt="%s")


db = DataBaseType.GSE52588
geos = [GeoType.any, GeoType.islands, GeoType.islands_shores]

for geo in geos:

    config = Config(
        db=db,
        geo=geo
    )

    save_gene_data(config)
