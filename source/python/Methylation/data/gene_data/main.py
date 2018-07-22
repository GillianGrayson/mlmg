from config.config import *
from config.types import *
from data.gene_data.gene_data import save_gene_data

config = Config(
    db=DataBaseType.GSE40279,
    geo=GeoType.islands
)

save_gene_data(config)
