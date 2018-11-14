from config.config import *
from config.types.annotations import GeneDataType, GeoType
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
import pandas as pd

num_top = 500

method = Method.manova

config_f = Config(
    read_only=True,
    db=DataBase.GSE40279,
    dt=DataType.bop,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=method,
    gender=Gender.F,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

config_m = Config(
    read_only=True,
    db=DataBase.GSE40279,
    dt=DataType.bop,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=method,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

keys = ['bop', 'gene', 'pval']

f_top_dict = load_top_dict(config_f, keys)
m_top_dict = load_top_dict(config_m, keys)

f_gene_names = list(f_top_dict['gene'])
m_gene_names = list(m_top_dict['gene'])

i_gene_names = list(set(f_gene_names).intersection(m_gene_names))
only_f_gene_names = list(set(f_gene_names) - set(i_gene_names))
only_m_gene_names = list(set(m_gene_names) - set(i_gene_names))

f_bop_names = list(f_top_dict['bop'])
m_bop_names = list(m_top_dict['bop'])

i_bop_names = list(set(f_bop_names).intersection(m_bop_names))
only_f_bop_names = list(set(f_bop_names) - set(i_bop_names))
only_m_bop_names = list(set(m_bop_names) - set(i_bop_names))

only_f_ids = []
only_f_bop = []
only_f_gene = []
only_f_p_value = []

for f_id in range(0, len(only_f_bop_names)):
    name = only_f_bop_names[f_id]
    index = f_bop_names.index(name)
    only_f_ids.append(index)
    only_f_gene.append(f_top_dict['gene'][index])
    only_f_p_value.append(f_top_dict['pval'][index])

only_f_order = np.argsort(only_f_ids)
only_f_bop_names_sorted = list(np.array(only_f_bop_names)[only_f_order])
only_f_ids_sorted = list(np.array(only_f_ids)[only_f_order])
only_f_gene_sorted = list(np.array(only_f_gene)[only_f_order])
only_f_p_value_sorted = list(np.array(only_f_p_value)[only_f_order])

strs = []
for f_id in range(0, len(only_f_bop_names)):
    name = only_f_bop_names_sorted[f_id]
    curr_str = name
    curr_str += (' ' + str(format(only_f_ids_sorted[f_id] + 1, 'd')))
    curr_str += (' ' + only_f_gene_sorted[f_id])
    curr_str += (' ' + str(format(float(only_f_p_value_sorted[f_id]), '0.4e')))
    strs.append(curr_str)
fn = 'only_f.txt'
#np.savetxt(fn, strs, fmt="%s")

only_m_ids = []
only_m_bop = []
only_m_gene = []
only_m_p_value = []

for m_id in range(0, len(only_m_bop_names)):
    name = only_m_bop_names[m_id]
    index = m_bop_names.index(name)
    only_m_ids.append(index)
    only_m_gene.append(m_top_dict['gene'][index])
    only_m_p_value.append(m_top_dict['pval'][index])

only_m_order = np.argsort(only_m_ids)
only_m_bop_names_sorted = list(np.array(only_m_bop_names)[only_m_order])
only_m_ids_sorted = list(np.array(only_m_ids)[only_m_order])
only_m_gene_sorted = list(np.array(only_m_gene)[only_m_order])
only_m_p_value_sorted = list(np.array(only_m_p_value)[only_m_order])

strs = []
for m_id in range(0, len(only_m_bop_names)):
    name = only_m_bop_names_sorted[m_id]
    curr_str = name
    curr_str += (' ' + str(format(only_m_ids_sorted[m_id] + 1, 'd')))
    curr_str += (' ' + only_m_gene_sorted[m_id])
    curr_str += (' ' + str(format(float(only_m_p_value_sorted[m_id]), '0.4e')))
    strs.append(curr_str)
fn = 'only_m.txt'
#np.savetxt(fn, strs, fmt="%s")

i_f_ids = []
i_m_ids = []
i_f_bop = []
i_m_bop = []
i_f_gene = []
i_m_gene = []
i_f_p_value = []
i_m_p_value = []

for i_id in range(0, len(i_bop_names)):
    name = i_bop_names[i_id]
    f_index = f_bop_names.index(name)
    m_index = m_bop_names.index(name)
    i_f_ids.append(f_index)
    i_m_ids.append(m_index)
    i_f_gene.append(f_top_dict['gene'][f_index])
    i_m_gene.append(m_top_dict['gene'][m_index])
    i_f_p_value.append(f_top_dict['pval'][f_index])
    i_m_p_value.append(m_top_dict['pval'][m_index])

i_order = np.argsort(i_f_ids)
i_bop_names_sorted = list(np.array(i_bop_names)[i_order])
i_f_ids_sorted = list(np.array(i_f_ids)[i_order])
i_m_ids_sorted = list(np.array(i_m_ids)[i_order])
i_f_gene_sorted = list(np.array(i_f_gene)[i_order])
i_m_gene_sorted = list(np.array(i_m_gene)[i_order])
i_f_p_value_sorted = list(np.array(i_f_p_value)[i_order])
i_m_p_value_sorted = list(np.array(i_m_p_value)[i_order])

strs = []
for i_id in range(0, len(i_bop_names)):
    name = i_bop_names_sorted[i_id]
    curr_str = name
    curr_str += (' ' + str(format(i_f_ids_sorted[i_id] + 1, 'd')))
    curr_str += (' ' + str(format(i_m_ids_sorted[i_id] + 1, 'd')))
    curr_str += (' ' + i_f_gene_sorted[i_id])
    curr_str += (' ' + i_m_gene_sorted[i_id])
    curr_str += (' ' + str(format(float(i_m_p_value_sorted[i_id]), '0.4e')))
    curr_str += (' ' + str(format(float(i_f_p_value_sorted[i_id]), '0.4e')))
    strs.append(curr_str)
fn = 'i.txt'
#np.savetxt(fn, strs, fmt="%s")

df_F = pd.DataFrame({'F Gene': only_f_gene_sorted,
                     'F Top': [x+1 for x in only_f_ids_sorted],
                     'F Bop': only_f_bop_names_sorted,
                     'F P_val': only_f_p_value_sorted})
df_I = pd.DataFrame({'I Gene': i_f_gene_sorted,
                     'I Top_F': [x+1 for x in i_f_ids_sorted],
                     'I Top_M': [x+1 for x in i_m_ids_sorted],
                     'I Bop': i_bop_names_sorted,
                     'I_P_val_F': i_f_p_value_sorted,
                     'I_P_val_M': i_m_p_value_sorted})
df_M = pd.DataFrame({'M Gene': only_m_gene_sorted,
                     'M Top': [x+1 for x in only_m_ids_sorted],
                     'M Bop': only_m_bop_names_sorted,
                     'M P_val': only_m_p_value_sorted})

writer = pd.ExcelWriter('pandas_simple.xlsx', engine='xlsxwriter')
df_F.to_excel(writer, index=False, sheet_name='F')
df_M.to_excel(writer, index=False, sheet_name='M')
df_I.to_excel(writer, index=False, sheet_name='I')
writer.save()
