from config.config import *
from config.types.annotations import GeneDataType, GeoType, ChromosomeType
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from config.types.experiments.method import *
from infrastructure.save.features import save_features
from gene.validation.top.gender_specific.butterfly.butterfly_dict import get_butterfly_dict


def data_bases_versus(config, data_bases, part):
    butterfly_dicts = []
    for data_base in data_bases:
        config = Config(
            read_only=True,
            data_base=data_base,
            data_type=config.data_type,

            chromosome_type=config.chromosome_type,

            geo_type=config.geo_type,
            gene_data_type=config.gene_data_type,

            disease=config.disease,

            scenario=config.scenario,
            approach=config.approach,
            method=config.method,

            is_clustering = config.is_clustering
        )

        butterfly_dicts.append(get_butterfly_dict(config))

    intersection_genes = butterfly_dicts[0]['gene']
    num_genes = int(len(intersection_genes) * part)
    for butterfly_dict in butterfly_dicts:
        genes = butterfly_dict['gene'][0:num_genes]
        intersection_genes = list(set(intersection_genes).intersection(genes))

    for gene in intersection_genes:
        print(gene)
    print(len(intersection_genes))

    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)

    fn = 'intersection_butterfly_genes_data_bases(' + data_bases_str + ')_method(' + config.method.value + ')_part(' + str(part) + ').txt'
    config_dump = Config(
        read_only=True,
        data_base=DataBase.data_base_versus,
        data_type=config.data_type,

        chromosome_type=config.chromosome_type,

        geo_type=config.geo_type,
        gene_data_type=config.gene_data_type,

        disease=config.disease,
        gender=Gender.versus,

        scenario=Scenario.validation,
        approach=Approach.top,
        method=Method.gender_specific,

        is_clustering = config.is_clustering
    )

    features = [
        intersection_genes
    ]

    fn = get_result_path(config_dump, fn)
    save_features(fn, features)


part = 0.05

data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_type = DataType.gene

chromosome_type = ChromosomeType.non_gender

gene_data_type = GeneDataType.mean
geo_type = GeoType.islands_shores

disease = Disease.any

scenario = Scenario.approach
approach = Approach.top
method = Method.linreg

is_clustering = False

config = Config(
    read_only=True,

    data_type=data_type,

    chromosome_type=chromosome_type,

    gene_data_type=gene_data_type,
    geo_type=geo_type,

    disease=disease,

    scenario=scenario,
    approach=approach,
    method=method,

    is_clustering=is_clustering
)

data_bases_versus(config, data_bases, part)
