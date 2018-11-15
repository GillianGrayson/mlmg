from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from config.types.experiments.method import *
from annotations.gene import get_dict_cpg_gene
from infrastructure.save.features import save_features
from cpg.validation.top.gender_specific.butterfly.butterfly_dict import get_butterfly_dict


def data_bases_versus(config, data_bases, part):
    butterfly_dicts = []
    for data_base in data_bases:
        config = Config(
            read_only=True,

            data_base=data_base,
            data_type=config.data_type,

            chromosome_type=config.chromosome_type,

            dna_region=config.dna_region,

            scenario=config.scenario,
            approach=config.approach,
            method=config.method,

            disease=config.disease,

            is_clustering = config.is_clustering
        )

        butterfly_dicts.append(get_butterfly_dict(config))

    intersection_cpgs = butterfly_dicts[0]['cpg']
    num_genes = int(len(intersection_cpgs) * part)
    for butterfly_dict in butterfly_dicts:
        cpgs = butterfly_dict['cpg'][0:num_genes]
        intersection_cpgs = list(set(intersection_cpgs).intersection(cpgs))

    dict_cpg_gene = get_dict_cpg_gene(config)

    intersection_genes = []
    for cpg in intersection_cpgs:
        print(cpg)
        if cpg in dict_cpg_gene:
            genes = dict_cpg_gene.get(cpg)
            for gene in genes:
                if gene not in intersection_genes:
                    intersection_genes.append(gene)
    print(len(intersection_cpgs))

    for gene in intersection_genes:
        print(gene)
    print(len(intersection_genes))

    data_bases_str = [x.value for x in data_bases]
    data_bases_str.sort()
    data_bases_str = '_'.join(data_bases_str)


    config_dump = Config(
        read_only=True,
        data_base=DataBase.data_base_versus,
        data_type=config.data_type,

        chromosome_type=config.chromosome_type,

        dna_region=config.dna_region,

        disease=config.disease,
        gender=Gender.versus,

        scenario=Scenario.validation,
        approach=Approach.top,
        method=Method.gender_specific,

        is_clustering = config.is_clustering
    )

    features = [
        intersection_cpgs
    ]
    fn = 'intersection_butterfly_cpgs_data_bases(' + data_bases_str + ')_method(' + config.method.value + ')_part(' + str(part) + ').txt'
    fn = get_result_path(config_dump, fn)
    save_features(fn, features)

    features = [
        intersection_genes
    ]
    fn = 'intersection_butterfly_genes_data_bases(' + data_bases_str + ')_method(' + config.method.value + ')_part(' + str(part) + ').txt'
    fn = get_result_path(config_dump, fn)
    save_features(fn, features)


part = 0.05

data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_type = DataType.cpg

chromosome_type = ChromosomeType.non_gender

dna_region = DNARegionType.genic

disease = Disease.any

scenario = Scenario.approach
approach = Approach.top
method = Method.linreg

is_clustering = False

config = Config(
    read_only=True,

    data_type=data_type,

    chromosome_type=chromosome_type,

    dna_region=dna_region,

    scenario=scenario,
    approach=approach,
    method=method,

    disease=disease,

    is_clustering=is_clustering
)

data_bases_versus(config, data_bases, part)
