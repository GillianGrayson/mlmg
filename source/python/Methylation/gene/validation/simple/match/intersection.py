from config.config import *
from infrastructure.load.top import *
import numpy as np

num_top = 500

# configs = [
#     Config(
#         db=DataBaseType.GSE40279,
#         dt=DataType.gene,
#         scenario=Scenario.approach,
#         approach=Approach.top,
#         approach_method=Method.manova,
#         gender=Gender.any,
#         disease=Disease.any,
#         approach_gd=GeneDataType.from_bop
#     ),
#     Config(
#         db=DataBaseType.GSE52588,
#         dt=DataType.gene,
#         scenario=Scenario.approach,
#         approach=Approach.top,
#         approach_method=Method.manova,
#         gender=Gender.any,
#         disease=Disease.healthy,
#         approach_gd=GeneDataType.from_bop
#     ),
#     Config(
#         db=DataBaseType.GSE30870,
#         dt=DataType.gene,
#         scenario=Scenario.approach,
#         approach=Approach.top,
#         approach_method=Method.manova,
#         gender=Gender.any,
#         disease=Disease.any,
#         approach_gd=GeneDataType.from_bop
#     )
# ]

configs = [
    Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.manova,
        gender=Gender.F,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    ),
    Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.linreg,
        gender=Gender.F,
        disease=Disease.any,
        approach_gd=GeneDataType.from_cpg
    ),
    Config(
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
]

intersection_genes = load_top_gene_names(configs[0], num_top)[0:num_top]
for config in configs[1:]:
    genes = load_top_gene_names(config, num_top)[0:num_top]
    intersection_genes = list(set(intersection_genes).intersection(genes))
print(len(intersection_genes))
genes_f = intersection_genes


configs = [
    Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.manova,
        gender=Gender.M,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    ),
    Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.linreg,
        gender=Gender.M,
        disease=Disease.any,
        approach_gd=GeneDataType.from_cpg
    ),
    Config(
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
]

intersection_genes = load_top_gene_names(configs[0], num_top)[0:num_top]
for config in configs[1:]:
    genes = load_top_gene_names(config, num_top)[0:num_top]
    intersection_genes = list(set(intersection_genes).intersection(genes))
print(len(intersection_genes))
genes_m = intersection_genes

intersection_genes = list(set(genes_m).intersection(genes_f))
only_m = list(set(genes_m) - set(intersection_genes))
only_f = list(set(genes_f) - set(intersection_genes))

fn = 'hannum2013_genes.txt'
article = load_top_gene_names_by_article(configs[0], fn)

intersection_genes.sort()
print('intersection genes:')
for gene in intersection_genes:
    print(gene)
print('\n\n')
intersection_article = list(set(intersection_genes).intersection(article))
intersection_article.sort()
for gene in intersection_article:
    print(gene)
print('\n\n')

only_m.sort()
print('only_m genes:')
for gene in only_m:
    print(gene)
print('\n\n')
intersection_article = list(set(only_m).intersection(article))
intersection_article.sort()
for gene in intersection_article:
    print(gene)
print('\n\n')

only_f.sort()
print('only_f genes:')
for gene in only_f:
    print(gene)
print('\n\n')
intersection_article = list(set(only_f).intersection(article))
intersection_article.sort()
for gene in intersection_article:
    print(gene)
print('\n\n')