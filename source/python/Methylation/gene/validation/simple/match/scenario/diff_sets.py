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
    approach_method=Method.linreg,
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
    approach_method=Method.linreg,
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

A_M_genes = load_top_gene_names(A_M, num_top)[0:num_top]
A_F_genes = load_top_gene_names(A_F, num_top)[0:num_top]
A_I_genes = list(set(A_M_genes).intersection(A_F_genes))
A_M_genes = list(set(A_M_genes) - set(A_I_genes))
A_F_genes = list(set(A_F_genes) - set(A_I_genes))

B_M_genes = load_top_gene_names(B_M, num_top)[0:num_top]
B_F_genes = load_top_gene_names(B_F, num_top)[0:num_top]
B_I_genes = list(set(B_M_genes).intersection(B_F_genes))
B_M_genes = list(set(B_M_genes) - set(B_I_genes))
B_F_genes = list(set(B_F_genes) - set(B_I_genes))

C_M_genes = load_top_gene_names(C_M, num_top)[0:num_top]
C_F_genes = load_top_gene_names(C_F, num_top)[0:num_top]
C_I_genes = list(set(C_M_genes).intersection(C_F_genes))
C_M_genes = list(set(C_M_genes) - set(C_I_genes))
C_F_genes = list(set(C_F_genes) - set(C_I_genes))

green_M_genes = A_M_genes
green_M_genes = list(set(green_M_genes).intersection(B_M_genes))
green_M_genes = list(set(green_M_genes).intersection(C_M_genes))

green_F_genes = A_F_genes
green_F_genes = list(set(green_F_genes).intersection(B_F_genes))
green_F_genes = list(set(green_F_genes).intersection(C_F_genes))

green_I_genes = A_I_genes
green_I_genes = list(set(green_I_genes).intersection(B_I_genes))
green_I_genes = list(set(green_I_genes).intersection(C_I_genes))

red_M_genes = B_M_genes
red_M_genes = list(set(red_M_genes).intersection(A_M_genes))
only_red_M_genes = list(set(red_M_genes) - set(green_M_genes))

red_F_genes = B_F_genes
red_F_genes = list(set(red_F_genes).intersection(A_F_genes))
only_red_F_genes = list(set(red_F_genes) - set(green_F_genes))

red_I_genes = B_I_genes
red_I_genes = list(set(red_I_genes).intersection(A_I_genes))
only_red_I_genes = list(set(red_I_genes) - set(green_I_genes))

blue_M_genes = B_M_genes
blue_M_genes = list(set(blue_M_genes).intersection(C_M_genes))
only_blue_M_genes = list(set(blue_M_genes) - set(green_M_genes))

blue_F_genes = B_F_genes
blue_F_genes = list(set(blue_F_genes).intersection(C_F_genes))
only_blue_F_genes = list(set(blue_F_genes) - set(green_F_genes))

blue_I_genes = B_I_genes
blue_I_genes = list(set(blue_I_genes).intersection(C_I_genes))
only_blue_I_genes = list(set(blue_I_genes) - set(green_I_genes))

only_black_M_genes = B_M_genes
only_black_M_genes = list(set(only_black_M_genes) - set(green_M_genes))
only_black_M_genes = list(set(only_black_M_genes) - set(only_red_M_genes))
only_black_M_genes = list(set(only_black_M_genes) - set(only_blue_M_genes))

only_black_F_genes = B_F_genes
only_black_F_genes = list(set(only_black_F_genes) - set(green_F_genes))
only_black_F_genes = list(set(only_black_F_genes) - set(only_red_F_genes))
only_black_F_genes = list(set(only_black_F_genes) - set(only_blue_F_genes))

only_black_I_genes = B_I_genes
only_black_I_genes = list(set(only_black_I_genes) - set(green_I_genes))
only_black_I_genes = list(set(only_black_I_genes) - set(only_red_I_genes))
only_black_I_genes = list(set(only_black_I_genes) - set(only_blue_I_genes))


green_M_genes.sort()
print('green_M_genes' + ' (' + str(len(green_M_genes)) + '):')
for gene in green_M_genes:
    print(gene)
print('\n\n')

only_red_M_genes.sort()
print('only_red_M_genes' + ' (' + str(len(only_red_M_genes)) + '):')
for gene in only_red_M_genes:
    print(gene)
print('\n\n')

only_blue_M_genes.sort()
print('only_blue_M_genes' + ' (' + str(len(only_blue_M_genes)) + '):')
for gene in only_blue_M_genes:
    print(gene)
print('\n\n')

only_black_M_genes.sort()
print('only_black_M_genes' + ' (' + str(len(only_black_M_genes)) + '):')
for gene in only_black_M_genes:
    print(gene)
print('\n\n')


green_F_genes.sort()
print('green_F_genes' + ' (' + str(len(green_F_genes)) + '):')
for gene in green_F_genes:
    print(gene)
print('\n\n')

only_red_F_genes.sort()
print('only_red_F_genes' + ' (' + str(len(only_red_F_genes)) + '):')
for gene in only_red_F_genes:
    print(gene)
print('\n\n')

only_blue_F_genes.sort()
print('only_blue_F_genes' + ' (' + str(len(only_blue_F_genes)) + '):')
for gene in only_blue_F_genes:
    print(gene)
print('\n\n')

only_black_F_genes.sort()
print('only_black_F_genes' + ' (' + str(len(only_black_F_genes)) + '):')
for gene in only_black_F_genes:
    print(gene)
print('\n\n')


green_I_genes.sort()
print('green_I_genes' + ' (' + str(len(green_I_genes)) + '):')
for gene in green_I_genes:
    print(gene)
print('\n\n')

only_red_I_genes.sort()
print('only_red_I_genes' + ' (' + str(len(only_red_I_genes)) + '):')
for gene in only_red_I_genes:
    print(gene)
print('\n\n')

only_blue_I_genes.sort()
print('only_blue_I_genes' + ' (' + str(len(only_blue_I_genes)) + '):')
for gene in only_blue_I_genes:
    print(gene)
print('\n\n')

only_black_I_genes.sort()
print('only_black_I_genes' + ' (' + str(len(only_black_I_genes)) + '):')
for gene in only_black_I_genes:
    print(gene)
