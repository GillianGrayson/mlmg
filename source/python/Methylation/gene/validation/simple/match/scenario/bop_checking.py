from config.config import *
from infrastructure.load.top import *

num_top = 500

D1 = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.bop,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.manova,
    gender=Gender.any,
    disease=Disease.any,
)

D2 = Config(
    read_only=True,
    db=DataBaseType.GSE30870,
    dt=DataType.bop,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.manova,
    gender=Gender.any,
    disease=Disease.any,
)

D3 = Config(
    read_only=True,
    db=DataBaseType.GSE52588,
    dt=DataType.bop,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=Method.manova,
    gender=Gender.any,
    disease=Disease.any,
)

D1_genes = load_top_data(D1, num_top, 1)
D2_genes = load_top_data(D2, num_top, 1)
D3_genes = load_top_data(D3, num_top, 1)

I_genes = list(set(list(set(D1_genes).intersection(D2_genes))).intersection(D3_genes))

fn = 'claudio2015_genes.txt'
article_genes = load_top_gene_names_by_article(Config(read_only=True), fn)

I_article_genes = list(set(I_genes).intersection(article_genes))


ololo = 1