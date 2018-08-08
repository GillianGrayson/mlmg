from config.config import *
from infrastructure.load.top import *


def vs_article(article_name, intersection, only_f, only_m):
    fn = article_name + '_genes.txt'
    article = load_top_gene_names_by_article(config_tmp, fn)
    intersection_article = list(set(intersection).intersection(article))
    intersection_article.sort()
    print('intersection_' + article_name + ' (' + str(len(intersection_article)) + '):')
    for gene in intersection_article:
       print(gene)
    print('\n\n')
    only_f_article = list(set(only_f).intersection(article))
    only_f_article.sort()
    print('only_f_' + article_name + ' (' + str(len(only_f_article)) + '):')
    for gene in only_f_article:
        print(gene)
    print('\n\n')
    only_m_article = list(set(only_m).intersection(article))
    only_m_article.sort()
    print('only_m_' + article_name + ' (' + str(len(only_m_article)) + '):')
    for gene in only_m_article:
        print(gene)
    print('\n\n')


num_top = 500

config_tmp = Config(read_only=True)

A_F = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.linreg,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.from_cpg
)

B_F = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.enet,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

C_F = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.manova,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.from_bop
)

D_F = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.random_forest,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

configs_F = [B_F]

I_F = load_top_gene_names(configs_F[0], num_top)[0:num_top]
for config in configs_F[1:]:
    genes = load_top_gene_names(config, num_top)[0:num_top]
    I_F = list(set(I_F).intersection(genes))
print('I_F len: ', len(I_F))

A_M = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.linreg,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.from_cpg
)

B_M = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.enet,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

C_M = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.manova,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.from_bop
)

D_M = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.random_forest,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

configs_M = [B_M]

I_M = load_top_gene_names(configs_M[0], num_top)[0:num_top]
for config in configs_M[1:]:
    genes = load_top_gene_names(config, num_top)[0:num_top]
    I_M = list(set(I_M).intersection(genes))
print('I_M len: ', len(I_M))

intersection = list(set(I_M).intersection(I_F))
intersection.sort()
print('intersection' + ' (' + str(len(intersection)) + '):')
for gene in intersection:
    print(gene)
print('\n\n')

only_f = list(set(I_F) - set(intersection))
only_f.sort()
print('only_f' + ' (' + str(len(only_f)) + '):')
for gene in only_f:
    print(gene)
print('\n\n')

only_m = list(set(I_M) - set(intersection))
only_m.sort()
print('only_m' + ' (' + str(len(only_m)) + '):')
for gene in only_m:
    print(gene)
print('\n\n')

article_name = 'claudio2015'
vs_article(article_name, intersection, only_f, only_m)

article_name = 'hannum2013'
vs_article(article_name, intersection, only_f, only_m)

article_name = 'singmann2015'
vs_article(article_name, intersection, only_f, only_m)

article_name = 'inoshita2015_cpg_cell_proportional'
vs_article(article_name, intersection, only_f, only_m)

article_name = 'inoshita2015_cpg_cell_nonproportional'
vs_article(article_name, intersection, only_f, only_m)

