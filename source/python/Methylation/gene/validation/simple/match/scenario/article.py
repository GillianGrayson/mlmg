from config.config import *
from config.types.annotations import GeneDataType, GeoType
from config.types.attributes.common import Disease, Gender
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
    db=DataBase.GSE40279,
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
    db=DataBase.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.linreg,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

C_F = Config(
    read_only=True,
    db=DataBase.GSE40279,
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
    db=DataBase.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.random_forest,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

A_M = Config(
    read_only=True,
    db=DataBase.GSE40279,
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
    db=DataBase.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.linreg,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

C_M = Config(
    read_only=True,
    db=DataBase.GSE40279,
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
    db=DataBase.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.random_forest,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)


A_M_genes = load_top_gene_names(A_M, num_top)[0:num_top]
B_M_genes = load_top_gene_names(B_M, num_top)[0:num_top]
C_M_genes = load_top_gene_names(C_M, num_top)[0:num_top]

A_F_genes = load_top_gene_names(A_F, num_top)[0:num_top]
B_F_genes = load_top_gene_names(B_F, num_top)[0:num_top]
C_F_genes = load_top_gene_names(C_F, num_top)[0:num_top]

M_genes = list(set(A_M_genes).intersection(B_M_genes).intersection(C_M_genes))
F_genes = list(set(A_F_genes).intersection(B_F_genes).intersection(C_F_genes))
I_genes = list(set(M_genes).intersection(F_genes))
M_genes = list(set(M_genes) - set(I_genes))
F_genes = list(set(F_genes) - set(I_genes))

article_name = 'singmann2015_genes.txt'
article = load_top_gene_names_by_article(config_tmp, article_name)
M_I_genes_1 = list(set(M_genes).intersection(article))
F_I_genes_1 = list(set(F_genes).intersection(article))
I_I_genes_1 = list(set(I_genes).intersection(article))
only_M_genes_1 = list(set(M_genes) - set(M_I_genes_1))
only_F_genes_1 = list(set(F_genes) - set(F_I_genes_1))
only_I_genes_1 = list(set(M_genes) - set(M_I_genes_1))

article_name = 'inoshita2015_cpg_cell_nonproportional_genes.txt'
article = load_top_gene_names_by_article(config_tmp, article_name)
M_I_genes_2 = list(set(M_genes).intersection(article))
F_I_genes_2 = list(set(F_genes).intersection(article))
I_I_genes_2 = list(set(I_genes).intersection(article))
only_M_genes_2 = list(set(M_genes) - set(M_I_genes_2))
only_F_genes_2 = list(set(F_genes) - set(F_I_genes_2))
only_I_genes_2 = list(set(M_genes) - set(M_I_genes_2))

M_green = list(set(M_I_genes_1).intersection(M_I_genes_2))
M_red =  list(set(M_I_genes_1) - set(M_green))
M_blue =  list(set(M_I_genes_2) - set(M_green))
M_black = list(set(M_genes) - set(M_green) - set(M_red) - set(M_blue))

M_green.sort()
print('M_green' + ' (' + str(len(M_green)) + '):')
for gene in M_green:
    print(gene)
print('\n\n')

M_red.sort()
print('M_red' + ' (' + str(len(M_red)) + '):')
for gene in M_red:
    print(gene)
print('\n\n')

M_blue.sort()
print('M_blue' + ' (' + str(len(M_blue)) + '):')
for gene in M_blue:
    print(gene)
print('\n\n')

M_black.sort()
print('M_black' + ' (' + str(len(M_black)) + '):')
for gene in M_black:
    print(gene)
print('\n\n')

F_green = list(set(F_I_genes_1).intersection(F_I_genes_2))
F_red =  list(set(F_I_genes_1) - set(F_green))
F_blue =  list(set(F_I_genes_2) - set(F_green))
F_black = list(set(F_genes) - set(F_green) - set(F_red) - set(F_blue))

F_green.sort()
print('F_green' + ' (' + str(len(F_green)) + '):')
for gene in F_green:
    print(gene)
print('\n\n')

F_red.sort()
print('F_red' + ' (' + str(len(F_red)) + '):')
for gene in F_red:
    print(gene)
print('\n\n')

F_blue.sort()
print('F_blue' + ' (' + str(len(F_blue)) + '):')
for gene in F_blue:
    print(gene)
print('\n\n')

F_black.sort()
print('F_black' + ' (' + str(len(F_black)) + '):')
for gene in F_black:
    print(gene)
print('\n\n')

I_green = list(set(I_I_genes_1).intersection(I_I_genes_2))
I_red =  list(set(I_I_genes_1) - set(I_green))
I_blue =  list(set(I_I_genes_2) - set(I_green))
I_black = list(set(I_genes) - set(I_green) - set(I_red) - set(I_blue))

I_green.sort()
print('I_green' + ' (' + str(len(I_green)) + '):')
for gene in I_green:
    print(gene)
print('\n\n')

I_red.sort()
print('I_red' + ' (' + str(len(I_red)) + '):')
for gene in I_red:
    print(gene)
print('\n\n')

I_blue.sort()
print('I_blue' + ' (' + str(len(I_blue)) + '):')
for gene in I_blue:
    print(gene)
print('\n\n')

I_black.sort()
print('I_black' + ' (' + str(len(I_black)) + '):')
for gene in I_black:
    print(gene)
print('\n\n')

article_name = 'inoshita2015_cpg_cell_proportional_genes.txt'
article = load_top_gene_names_by_article(config_tmp, article_name)

article_name = 'inoshita2015_cpg_cell_nonproportional_genes.txt'
article = load_top_gene_names_by_article(config_tmp, article_name)