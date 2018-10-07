from config.config import *
from infrastructure.load.top import *
from config.method import *
import pandas as pd
import math


rare_part = 0.005

data_base = DataBase.GSE87571
data_type = DataType.gene

chromosome_type = ChromosomeTypes.non_gender

gene_data_type = GeneDataType.mean
geo_type = GeoType.islands_shores

disease = Disease.any
scenario = Scenario.approach
approach = Approach.top
method = Method.linreg

is_clustering = False

method_metrics = get_method_metrics(method, is_clustering)

config_f = Config(
    read_only=True,
    data_base=data_base,
    data_type=data_type,

    chromosome_type=chromosome_type,

    geo_type=geo_type,
    gene_data_type=gene_data_type,

    disease=disease,
    gender=Gender.F,
    scenario=scenario,
    approach=approach,
    method=method
)

config_m = Config(
    read_only=True,
    data_base=data_base,
    data_type=data_type,

    chromosome_type=chromosome_type,

    geo_type=geo_type,
    gene_data_type=gene_data_type,

    disease=disease,
    gender=Gender.M,
    scenario=scenario,
    approach=approach,
    method=method
)

keys = ['gene'] + method_metrics

f_all_dict = load_top_dict(config_f, keys)
m_all_dict = load_top_dict(config_m, keys)

f_genes = f_all_dict['gene']
m_genes = m_all_dict['gene']

num_genes = len(f_genes)

f_metrics = f_all_dict[get_method_main_metric(method)]
m_metrics = m_all_dict[get_method_main_metric(method)]
for gene_id in range(0, num_genes):
    f_metrics[gene_id] = metric_processing(method, f_metrics[gene_id])
    m_metrics[gene_id] = metric_processing(method, m_metrics[gene_id])

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

num_rare = math.floor(rare_part * num_genes)

rare_genes = genes_srt[0:num_rare]
rare_diff = diff_srt[0:num_rare]

rare_dict = {'rare_genes':rare_genes, 'rare_diff': rare_diff}
rare_df = pd.DataFrame(rare_dict)

fn = 'rare_genes_' + 'method(' + method.value + ')_rare(' + str(rare_part) + ').xlsx'

config_f.scenario = Scenario.validation
config_f.method = Method.gender_specific
fn_f = get_result_path(config_f, fn)
writer = pd.ExcelWriter(fn_f, engine='xlsxwriter')
rare_df.to_excel(writer, index=False, sheet_name='rare_df')
writer.save()

config_m.scenario = Scenario.validation
config_m.method = Method.gender_specific
fn_m = get_result_path(config_m, fn)
writer = pd.ExcelWriter(fn_m, engine='xlsxwriter')
rare_df.to_excel(writer, index=False, sheet_name='rare_df')
writer.save()
