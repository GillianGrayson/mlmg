import numpy as np

from config import Config, ConfigGSE52588, ConfigGSE40279
from dicts import *
from geo import *

print_rate = 10000

fs_type = FSType.local_big
db_type = DataBaseType.GSE40279
geo_type = GeoType.any
config = Config(fs_type, db_type)
if db_type is DataBaseType.GSE40279:
    config = ConfigGSE40279(fs_type, db_type)
elif db_type is DataBaseType.GSE52588:
    config = ConfigGSE52588(fs_type, db_type)

gene_raw_dict = config.get_raw_dict(fs_type, db_type, geo_type)

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
    if num_genes % print_rate == 0:
        print('num_genes: ' + str(num_genes))

np.savetxt('gene_mean' + geo_type.value + '.txt', gene_mean_str_list, fmt="%s")
np.savetxt('gene_std' + geo_type.value + '.txt', gene_std_str_list, fmt="%s")
np.savetxt('gene_mean_der' + geo_type.value + '.txt', gene_mean_der_str_list, fmt="%s")
np.savetxt('gene_mean_der_normed' + geo_type.value + '.txt', gene_mean_der_normed_str_list, fmt="%s")


