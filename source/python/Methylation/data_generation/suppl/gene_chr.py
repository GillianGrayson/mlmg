from config.config import *
from annotations.gene import *


def save_gene_chr(config):
    gene_chr = get_dict_gene_chr(config)

    str_list = []
    for gene in gene_chr:
        curr_str = gene
        curr_gene_chr = gene_chr[gene]
        for chr in curr_gene_chr:
            curr_str += ' ' + str(chr)
        str_list.append(curr_str)

    fn = 'gene_chr.txt'
    fn = get_suppl_path(config, fn)
    np.savetxt(fn, str_list, fmt="%s")



config = Config(
    db=DataBase.GSE40279,
    geo=GeoType.any
)

save_gene_chr(config)
