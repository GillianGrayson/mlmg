from config.config import *
from infrastructure.load.top import *
from config.method import *
import pandas as pd
from infrastructure.save.features import save_features
import math


def get_butterfly_dict(config):

    method_metrics = get_method_metrics(config.method, config.is_clustering)

    config_f = Config(
        read_only=True,
        data_base=config.data_base,
        data_type=config.data_type,

        chromosome_type=config.chromosome_type,

        geo_type=config.geo_type,
        gene_data_type=config.gene_data_type,

        disease=config.disease,
        gender=Gender.F,
        scenario=config.scenario,
        approach=config.approach,
        method=config.method,

        is_clustering=config.is_clustering
    )

    config_m = Config(
        read_only=True,
        data_base=config.data_base,
        data_type=config.data_type,

        chromosome_type=config.chromosome_type,

        geo_type=config.geo_type,
        gene_data_type=config.gene_data_type,

        disease=config.disease,
        gender=Gender.M,
        scenario=config.scenario,
        approach=config.approach,
        method=config.method,

        is_clustering=config.is_clustering
    )

    keys = ['gene'] + method_metrics

    f_all_dict = load_top_dict(config_f, keys)
    m_all_dict = load_top_dict(config_m, keys)

    f_genes = f_all_dict['gene']
    m_genes = m_all_dict['gene']

    num_genes = len(f_genes)

    f_metrics = f_all_dict[get_method_main_metric(config.method)]
    m_metrics = m_all_dict[get_method_main_metric(config.method)]
    for gene_id in range(0, num_genes):
        f_metrics[gene_id] = metric_processing(config.method, f_metrics[gene_id])
        m_metrics[gene_id] = metric_processing(config.method, m_metrics[gene_id])

    f_metrics_passed = f_metrics
    m_metrics_passed = np.zeros(num_genes)
    metrics_diff = np.zeros(num_genes)
    for gene_id in range(0, num_genes):
        gene = f_genes[gene_id]
        m_id = m_genes.index(gene)
        m_metrics_passed[gene_id] = m_metrics[m_id]
        metrics_diff[gene_id] = f_metrics_passed[gene_id] - m_metrics_passed[gene_id]

    order = np.argsort(list(map(abs, metrics_diff)))[::-1]
    diff_srt = list(np.array(metrics_diff)[order])
    genes_srt = list(np.array(f_genes)[order])

    butterfly_dict = {'gene':genes_srt, 'diff': diff_srt}
    butterfly_df = pd.DataFrame(butterfly_dict)

    fn = 'butterfly_genes_' + 'method(' + config.method.value + ').xlsx'

    config_f.scenario = Scenario.validation
    config_f.method = Method.gender_specific
    fn_f = get_result_path(config_f, fn)
    writer = pd.ExcelWriter(fn_f, engine='xlsxwriter')
    butterfly_df.to_excel(writer, index=False, sheet_name='butterfly')
    writer.save()

    config_m.scenario = Scenario.validation
    config_m.method = Method.gender_specific
    fn_m = get_result_path(config_m, fn)
    writer = pd.ExcelWriter(fn_m, engine='xlsxwriter')
    butterfly_df.to_excel(writer, index=False, sheet_name='butterfly')
    writer.save()

    return butterfly_dict


def data_bases_intersection(config, data_bases, part):
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

    fn = 'intersection_genes_data_bases(' + data_bases_str + ')_method(' + config.method.value + ')_part(' + str(part) + ').txt'
    config_dump = Config(
        read_only=True,
        data_base=DataBase.data_base_versus,
        data_type=config.data_type,

        chromosome_type=config.chromosome_type,

        geo_type=config.geo_type,
        gene_data_type=config.gene_data_type,

        disease=config.disease,
        gender=Gender.F,

        scenario=Scenario.validation,
        approach=Approach.top,
        method=Method.gender_specific,

        is_clustering = config.is_clustering
    )

    features = [
        intersection_genes
    ]

    fn_f = get_result_path(config_dump, fn)
    save_features(fn_f, features)

    config_dump.gender = Gender.M
    fn_m = get_result_path(config_dump, fn)
    save_features(fn_m, features)


part = 0.05

data_bases = [DataBase.GSE87571, DataBase.GSE40279]
data_type = DataType.gene

chromosome_type = ChromosomeTypes.non_gender

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

data_bases_intersection(config, data_bases, part)
