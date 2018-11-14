from config.config import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from config.types.experiments.method import *
import pandas as pd

num_top = 500

data_base = DataBase.GSE87571
data_type = DataType.gene

chromosome_type = ChromosomeType.non_gender

gene_data_type = GeneDataType.mean
geo_type = GeoType.islands_shores

disease = Disease.any
scenario = Scenario.approach
approach = Approach.top
method = Method.linreg

method_metrics = get_method_metrics(method)

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

f_top_dict = load_top_dict(config_f, keys, num_top)
f_all_dict = load_top_dict(config_f, keys)
m_top_dict = load_top_dict(config_m, keys, num_top)
m_all_dict = load_top_dict(config_m, keys)

f_genes = f_top_dict['gene']
m_genes = m_top_dict['gene']

i_genes = list(set(f_genes).intersection(m_genes))
only_f_genes = list(set(f_genes) - set(i_genes))
only_m_genes = list(set(m_genes) - set(i_genes))

# Table for female-only
only_f_dict = {}
only_f_dict['gene'] = []
only_f_dict['top'] = []
only_f_dict['top_vs'] = []
for x in method_metrics:
    only_f_dict[x] = []
    only_f_dict[x + '_vs'] = []

for f_id in range(0, len(only_f_genes)):
    gene = only_f_genes[f_id]
    gene_id = f_top_dict['gene'].index(gene)
    is_in_vs = True if gene in m_all_dict['gene'] else False
    if is_in_vs:
        gene_id_vs = m_all_dict['gene'].index(gene)
    else:
        gene_id_vs = -1

    only_f_dict['gene'].append(gene)
    only_f_dict['top'].append(gene_id)
    if is_in_vs:
        only_f_dict['top_vs'].append(gene_id_vs)
    else:
        only_f_dict['top_vs'].append('None')
    for x in method_metrics:
        only_f_dict[x].append(f_top_dict[x][gene_id])
        if is_in_vs:
            only_f_dict[x + '_vs'].append(m_all_dict[x][gene_id_vs])
        else:
            only_f_dict[x + '_vs'].append('None')

only_f_order = np.argsort(only_f_dict['top'])
only_f_dict['gene'] = list(np.array(only_f_dict['gene'])[only_f_order])
only_f_dict['top'] = list(np.array(only_f_dict['top'])[only_f_order])
only_f_dict['top_vs'] = list(np.array(only_f_dict['top_vs'])[only_f_order])

for x in method_metrics:
    only_f_dict[x] = list(np.array(only_f_dict[x])[only_f_order])
    only_f_dict[x + '_vs'] = list(np.array(only_f_dict[x + '_vs'])[only_f_order])

only_f_df = pd.DataFrame(only_f_dict)

# Table for male-only
only_m_dict = {}
only_m_dict['gene'] = []
only_m_dict['top'] = []
only_m_dict['top_vs'] = []
for x in method_metrics:
    only_m_dict[x] = []
    only_m_dict[x + '_vs'] = []

for f_id in range(0, len(only_m_genes)):
    gene = only_m_genes[f_id]
    gene_id = m_top_dict['gene'].index(gene)
    is_in_vs = True if gene in f_all_dict['gene'] else False
    if is_in_vs:
        gene_id_vs = f_all_dict['gene'].index(gene)
    else:
        gene_id_vs = -1

    only_m_dict['gene'].append(gene)
    only_m_dict['top'].append(gene_id)
    if is_in_vs:
        only_m_dict['top_vs'].append(gene_id_vs)
    else:
        only_m_dict['top_vs'].append('None')
    for x in method_metrics:
        only_m_dict[x].append(m_top_dict[x][gene_id])
        if is_in_vs:
            only_m_dict[x + '_vs'].append(f_all_dict[x][gene_id_vs])
        else:
            only_m_dict[x + '_vs'].append('None')

only_m_order = np.argsort(only_m_dict['top'])
only_m_dict['gene'] = list(np.array(only_m_dict['gene'])[only_m_order])
only_m_dict['top'] = list(np.array(only_m_dict['top'])[only_m_order])
only_m_dict['top_vs'] = list(np.array(only_m_dict['top_vs'])[only_m_order])

for x in method_metrics:
    only_m_dict[x] = list(np.array(only_m_dict[x])[only_m_order])
    only_m_dict[x + '_vs'] = list(np.array(only_m_dict[x + '_vs'])[only_m_order])

only_m_df = pd.DataFrame(only_m_dict)

# Table for intersection
i_dict = {}
i_dict['gene'] = []
i_dict['top_f'] = []
i_dict['top_m'] = []
for x in method_metrics:
    i_dict[x + '_f'] = []
    i_dict[x + '_m'] = []

for f_id in range(0, len(i_genes)):
    gene = i_genes[f_id]
    gene_id_f = f_top_dict['gene'].index(gene)
    gene_id_m = m_top_dict['gene'].index(gene)

    i_dict['gene'].append(gene)
    i_dict['top_f'].append(gene_id_f)
    i_dict['top_m'].append(gene_id_m)
    for x in method_metrics:
        i_dict[x + '_f'].append(f_top_dict[x][gene_id_f])
        i_dict[x + '_m'].append(m_top_dict[x][gene_id_m])

i_order = np.argsort(i_dict['top_f'])
i_dict['gene'] = list(np.array(i_dict['gene'])[i_order])
i_dict['top_f'] = list(np.array(i_dict['top_f'])[i_order])
i_dict['top_m'] = list(np.array(i_dict['top_m'])[i_order])

for x in method_metrics:
    i_dict[x + '_f'] = list(np.array(i_dict[x + '_f'])[i_order])
    i_dict[x + '_m'] = list(np.array(i_dict[x + '_m'])[i_order])

i_df = pd.DataFrame(i_dict)

fn = 'sets_' + 'method(' + method.value + ')_top(' + str(num_top) + ').xlsx'
config_f.scenario = Scenario.validation
config_f.method = Method.gender_specific
fn_f = get_result_path(config_f, fn)
writer = pd.ExcelWriter(fn_f, engine='xlsxwriter')
only_f_df.to_excel(writer, index=False, sheet_name='only_f')
only_m_df.to_excel(writer, index=False, sheet_name='only_m')
i_df.to_excel(writer, index=False, sheet_name='i')
writer.save()

config_m.scenario = Scenario.validation
config_m.method = Method.gender_specific
fn_m = get_result_path(config_m, fn)
writer = pd.ExcelWriter(fn_m, engine='xlsxwriter')
only_f_df.to_excel(writer, index=False, sheet_name='only_f')
only_m_df.to_excel(writer, index=False, sheet_name='only_m')
i_df.to_excel(writer, index=False, sheet_name='i')
writer.save()

