from config.config import *
from infrastructure.load.top import *
from gene.validation.simple.match.top_dicts import get_config_dict
from gene.validation.simple.match.types import *
from gene.validation.simple.match.intersection import get_gene_intersection, get_gene_intersection_with_article


def gender_and_articles(fn, genes_I, gene_only_F, gene_only_M):
    article_intersection_genes_only_F = get_gene_intersection_with_article(fn, gene_only_F)
    print('genes_only_F with ' + fn + ' (' + str(len(article_intersection_genes_only_F)) + '):')
    for gene in article_intersection_genes_only_F:
        print(gene)
    print('\n\n')
    article_intersection_genes_only_M = get_gene_intersection_with_article(fn, gene_only_M)
    print('genes_only_M with ' + fn + ' (' + str(len(article_intersection_genes_only_M)) + '):')
    for gene in article_intersection_genes_only_M:
        print(gene)
    print('\n\n')
    article_intersection_genes_I = get_gene_intersection_with_article(fn, genes_I)
    print('genes_I with ' + fn + ' (' + str(len(article_intersection_genes_I)) + '):')
    for gene in article_intersection_genes_I:
        print(gene)
    print('\n\n')


num_top = 500

config_tmp = Config(read_only=True)
config_dict = get_config_dict(cpg_method=Method.enet,
                              gene_method=Method.enet,
                              bop_method=Method.manova)

configs_F = [config_dict[TopSource.GENE_F.value]]
genes_F = get_gene_intersection(configs_F, num_top)
print('genes_F len: ', len(genes_F))

configs_M = [config_dict[TopSource.GENE_M.value]]
genes_M = get_gene_intersection(configs_M, num_top)
print('genes_M len: ', len(genes_M))

genes_I = list(set(genes_F).intersection(genes_M))
genes_I_order = [genes_F.index(x) for x in genes_I]
genes_I_order.sort()
genes_I = list(np.array(genes_F)[genes_I_order])
print('genes_I' + ' (' + str(len(genes_I)) + '):')
for gene in genes_I:
    print(gene)
print('\n\n')

gene_only_F = list(set(genes_F) - set(genes_I))
gene_only_F_order = [genes_F.index(x) for x in gene_only_F]
gene_only_F_order.sort()
gene_only_F = list(np.array(genes_F)[gene_only_F_order])
print('gene_only_F' + ' (' + str(len(gene_only_F)) + '):')
for gene in gene_only_F:
    print(gene)
print('\n\n')

gene_only_M = list(set(genes_M) - set(genes_I))
gene_only_M_order = [genes_M.index(x) for x in gene_only_M]
gene_only_M_order.sort()
gene_only_M = list(np.array(genes_M)[gene_only_M_order])
print('gene_only_M' + ' (' + str(len(gene_only_M)) + '):')
for gene in gene_only_M:
    print(gene)
print('\n\n')

fn = 'claudio2015_genes.txt'
gender_and_articles(fn, genes_I, gene_only_F, gene_only_M)

fn = 'hannum2013_genes.txt'
gender_and_articles(fn, genes_I, gene_only_F, gene_only_M)

fn = 'singmann2015_genes.txt'
gender_and_articles(fn, genes_I, gene_only_F, gene_only_M)

article_name = 'inoshita2015_cpg_cell_proportional_genes.txt'
gender_and_articles(fn, genes_I, gene_only_F, gene_only_M)

fn = 'inoshita2015_cpg_cell_nonproportional_genes.txt'
gender_and_articles(fn, genes_I, gene_only_F, gene_only_M)

