from config.config import *
from infrastructure.load.top import *
from infrastructure.save.features import save_features


config = Config(
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    approach=Approach.top,
    scenario=Scenario.approach,
    approach_method=Method.enet,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands,
)

num_top = 100
fn = 'claudio2015.txt'

origin_top = load_top_gene_names_by_article(config, fn)
gene_names = load_top_gene_names(config, num_top)
genes_match = []
for gene in gene_names[0:num_top]:
    if gene in origin_top:
        genes_match.append(gene)

fn = 'match.txt'
fn = get_result_path(config, fn)
save_features(fn, genes_match)

print('top: ' + str(len(genes_match)))
