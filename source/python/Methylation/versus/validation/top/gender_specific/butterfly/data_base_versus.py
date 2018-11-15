from config.config import *
from config.types.annotations import GeneDataType, GeoType, ClassType, DNARegion, ChromosomeType
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from config.types.experiments.method import *


def data_types_genes_intersection(config, data_bases, data_types, part):
    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)
    gene_lists = []
    for data_type in data_types:
        method = data_types.get(data_type)
        curr_config = Config(
            read_only=True,

            data_base=config.data_base,
            data_type=data_type,

            chromosome_type=config.chromosome_type,

            dna_region=config.dna_region,

            class_type=config.class_type,

            geo_type=config.geo_type,
            gene_data_type=config.gene_data_type,

            scenario=config.scenario,
            approach=config.approach,
            method=config.method,

            disease=config.disease,
            gender=config.gender,

            is_clustering = config.is_clustering
        )

        fn = 'intersection_butterfly_genes_data_bases(' + data_bases_str + ')_method(' + method.value + ')_part(' + str(part) + ').txt'
        keys = ['gene']
        gene_dict = load_top_dict(curr_config, keys, fn=fn)
        gene_lists.append(gene_dict['gene'])


    gene_intersection = gene_lists[0]
    for genes in gene_lists:
        gene_intersection = list(set(gene_intersection).intersection(genes))

    for gene in gene_intersection:
        print(gene)
    print(str(len(gene_intersection)))


part = 0.05

data_base = DataBase.data_base_versus
data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_types = {DataType.gene: Method.linreg, DataType.bop: Method.manova, DataType.cpg:Method.linreg}

chromosome_type = ChromosomeType.non_gender

# bop
class_type = ClassType.class_ab
# cpg
dna_region = DNARegion.genic
# gene
geo_type = GeoType.islands_shores
gene_data_type = GeneDataType.mean

scenario = Scenario.validation
approach = Approach.top
method = Method.gender_specific

disease = Disease.any
gender = Gender.versus

is_clustering = False

config = Config(
    read_only=True,

    data_base=data_base,

    chromosome_type=chromosome_type,

    class_type=class_type,

    dna_region=dna_region,

    geo_type=geo_type,
    gene_data_type=gene_data_type,

    scenario=scenario,
    approach=approach,
    method=method,

    disease=disease,
    gender=gender,

    is_clustering=is_clustering
)

data_types_genes_intersection(config, data_bases, data_types, part)
