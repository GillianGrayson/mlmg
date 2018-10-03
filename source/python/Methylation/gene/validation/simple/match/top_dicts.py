from config.config import *
from infrastructure.load.top import *
from gene.validation.simple.match.types import *


def get_config_dict(cpg_method=Method.enet,
                    gene_method=Method.enet,
                    bop_method=Method.manova):

    config_dict = {}

    CPG_F = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=cpg_method,
        gender=Gender.F,
        disease=Disease.any,
        approach_gd=GeneDataType.from_cpg
    )

    GENE_F = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=gene_method,
        gender=Gender.F,
        disease=Disease.any,
        approach_gd=GeneDataType.mean,
        geo=GeoType.islands_shores
    )

    BOP_F = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=bop_method,
        gender=Gender.F,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    )

    CPG_M = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=cpg_method,
        gender=Gender.M,
        disease=Disease.any,
        approach_gd=GeneDataType.from_cpg
    )

    GENE_M = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=gene_method,
        gender=Gender.M,
        disease=Disease.any,
        approach_gd=GeneDataType.mean,
        geo=GeoType.islands_shores
    )

    BOP_M = Config(
        read_only=True,
        db=DataBase.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.approach,
        approach=Approach.top,
        approach_method=bop_method,
        gender=Gender.M,
        disease=Disease.any,
        approach_gd=GeneDataType.from_bop
    )

    config_dict[TopSource.CPG_F.value] = CPG_F
    config_dict[TopSource.GENE_F.value] = GENE_F
    config_dict[TopSource.BOP_F.value] = BOP_F

    config_dict[TopSource.CPG_M.value] = CPG_M
    config_dict[TopSource.GENE_M.value] = GENE_M
    config_dict[TopSource.BOP_M.value] = BOP_M

    return config_dict