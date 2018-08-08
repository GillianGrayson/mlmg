from config.config import *
from annotations.regular import *
from infrastructure.path import *
from infrastructure.save.features import save_features

def save_gene_by_cpg(config, fn):
    cpg_gene_dict = get_dict_cpg_gene(config)

    f = open(fn + '.txt')
    target_cpgs = f.read().splitlines()

    genes = []
    for cpg in target_cpgs:
        if cpg in cpg_gene_dict:
            curr_genes = cpg_gene_dict[cpg]
            for gene in curr_genes:
                if gene not in genes:
                    genes.append(gene)

    save_features(fn + '_genes.txt', genes)


config = Config(
    db=DataBaseType.GSE40279,
    geo=GeoType.any
)

fn = 'hannum2013_cpgs'
fn = get_path(config, fn)
save_gene_by_cpg(config, fn)
