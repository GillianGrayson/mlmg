from config.config import *
from config.types import *
from data.aux_data.aux_data import save_cpg_by_gene_list

config = Config(
    db=DataBaseType.GSE40279,
    geo=GeoType.any
)

fn = 'graph_genes_1'
save_cpg_by_gene_list(config, fn)
