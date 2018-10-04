import os
from infrastructure.directories.experiment import *
from infrastructure.tree.bop import create_bop_data_tree
from infrastructure.tree.gene import create_gene_data_tree
from infrastructure.tree.cpg import create_cpg_data_tree
from infrastructure.path.path import *


config = Config(data_base=DataBase.GSE87571,
                read_only=True)

create_bop_data_tree(config)
create_gene_data_tree(config)
create_cpg_data_tree(config)
