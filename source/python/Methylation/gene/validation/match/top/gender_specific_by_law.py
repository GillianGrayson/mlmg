from config.config import *
from infrastructure.load.top import *
from config.method import *
import pandas as pd
import math

num_top = 2000

data_base_type = DataBaseType.GSE87571
data_type = DataType.gene
scenario = Scenario.approach
approach = Approach.top

a_param = 58.12
b_param = 5.0

approach_method = Method.linreg
approach_method_metrics = get_method_metrics(approach_method)

disease = Disease.any
approach_gene_data_type = GeneDataType.mean
geo_type = GeoType.islands_shores

config_f = Config(
    read_only=True,
    db=data_base_type,
    dt=data_type,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=Gender.F,
    disease=disease,
    approach_gd=approach_gene_data_type,
    geo=geo_type
)

config_m = Config(
    read_only=True,
    db=data_base_type,
    dt=data_type,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=Gender.M,
    disease=disease,
    approach_gd=approach_gene_data_type,
    geo=geo_type
)

keys = ['gene'] + approach_method_metrics

f_all_dict = load_top_dict(config_f, keys)
m_all_dict = load_top_dict(config_m, keys)

f_genes = f_all_dict['gene']
m_genes = m_all_dict['gene']

# Table for female
f_genes_by_law = []
for id in range(0, num_top):
    gene = f_genes[id]
    id_vs = m_genes.index(gene)

    if approach_method is Method.linreg:
        allowed_distance = int(a_param * math.log(id + 1) + b_param)
    else:
        allowed_distance = 0

    if id_vs - id > allowed_distance:
        f_genes_by_law.append(gene)

f_dict = {}
f_dict['gene'] = []
f_dict['top'] = []
f_dict['top_vs'] = []
for x in approach_method_metrics:
    f_dict[x] = []
    f_dict[x + '_vs'] = []

for f_id in range(0, len(f_genes_by_law)):
    gene = f_genes_by_law[f_id]
    gene_id = f_all_dict['gene'].index(gene)
    is_in_vs = True if gene in m_all_dict['gene'] else False
    if is_in_vs:
        gene_id_vs = m_all_dict['gene'].index(gene)
    else:
        gene_id_vs = -1

    f_dict['gene'].append(gene)
    f_dict['top'].append(gene_id)
    if is_in_vs:
        f_dict['top_vs'].append(gene_id_vs)
    else:
        f_dict['top_vs'].append('None')
    for x in approach_method_metrics:
        f_dict[x].append(f_all_dict[x][gene_id])
        if is_in_vs:
            f_dict[x + '_vs'].append(m_all_dict[x][gene_id_vs])
        else:
            f_dict[x + '_vs'].append('None')

only_f_df = pd.DataFrame(f_dict)


# Table for male
m_genes_by_law = []
for id in range(0, num_top):
    gene = m_genes[id]
    id_vs = f_genes.index(gene)

    if approach_method is Method.linreg:
        allowed_distance = int(a_param * math.log(id + 1) + b_param)
    else:
        allowed_distance = 0

    if id_vs - id > allowed_distance:
        m_genes_by_law.append(gene)

m_dict = {}
m_dict['gene'] = []
m_dict['top'] = []
m_dict['top_vs'] = []
for x in approach_method_metrics:
    m_dict[x] = []
    m_dict[x + '_vs'] = []

for m_id in range(0, len(m_genes_by_law)):
    gene = m_genes_by_law[m_id]
    gene_id = m_all_dict['gene'].index(gene)
    is_in_vs = True if gene in f_all_dict['gene'] else False
    if is_in_vs:
        gene_id_vs = f_all_dict['gene'].index(gene)
    else:
        gene_id_vs = -1

    m_dict['gene'].append(gene)
    m_dict['top'].append(gene_id)
    if is_in_vs:
        m_dict['top_vs'].append(gene_id_vs)
    else:
        m_dict['top_vs'].append('None')
    for x in approach_method_metrics:
        m_dict[x].append(m_all_dict[x][gene_id])
        if is_in_vs:
            m_dict[x + '_vs'].append(f_all_dict[x][gene_id_vs])
        else:
            m_dict[x + '_vs'].append('None')

only_m_df = pd.DataFrame(m_dict)

writer = pd.ExcelWriter(data_base_type.value + '_' + approach_method.value +'_by_law_top(' + str(num_top) + ').xlsx', engine='xlsxwriter')
only_f_df.to_excel(writer, index=False, sheet_name='f')
only_m_df.to_excel(writer, index=False, sheet_name='m')
writer.save()
