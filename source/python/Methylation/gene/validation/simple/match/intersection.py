from config.config import *
from infrastructure.load.top import *
import numpy as np

num_top = 1000

configs = [
    Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.manova,
        gender=Gender.any,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    ),
    Config(
        db=DataBaseType.GSE52588,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.manova,
        gender=Gender.any,
        disease=Disease.healthy,
        approach_gd=GeneDataType.from_bop
    ),
    Config(
        db=DataBaseType.GSE30870,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=Method.manova,
        gender=Gender.any,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    )
]

intersection_genes = load_top_gene_names(configs[0], num_top)[0:num_top]

for config in configs[1:]:
    genes = load_top_gene_names(config, num_top)[0:num_top]
    intersection_genes = list(set(intersection_genes).intersection(genes))

fn = 'claudio2015_genes.txt'
genes = load_top_gene_names_by_article(configs[0], fn)
intersection_genes = list(set(intersection_genes).intersection(genes))

print(len(intersection_genes))

intersection_genes.sort()
for gene in intersection_genes:
    print(gene)

