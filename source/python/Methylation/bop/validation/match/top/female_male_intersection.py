from config.config import *
from infrastructure.load.top import *
from config.method import *
import pandas as pd

num_top = 500

data_base_type = DataBaseType.GSE40279
data_type = DataType.bop
scenario = Scenario.approach
approach = Approach.top

approach_method = Method.manova
approach_method_metrics = get_method_metrics(approach_method)

disease = Disease.any

config_f = Config(
    read_only=True,
    db=data_base_type,
    dt=data_type,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=Gender.F,
    disease=disease
)

config_m = Config(
    read_only=True,
    db=data_base_type,
    dt=data_type,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=Gender.M,
    disease=disease
)

keys = ['bop', 'gene'] + approach_method_metrics

f_top_dict = load_top_dict(config_f, keys, num_top)
f_all_dict = load_top_dict(config_f, keys)
m_top_dict = load_top_dict(config_m, keys, num_top)
m_all_dict = load_top_dict(config_m, keys)

f_bops = f_top_dict['bop']
m_bops = m_top_dict['bop']

i_bops = list(set(f_bops).intersection(m_bops))
only_f_bops = list(set(f_bops) - set(i_bops))
only_m_bops = list(set(m_bops) - set(i_bops))

# Table for female-only
only_f_dict = {}
only_f_dict['bop'] = []
only_f_dict['gene'] = []
only_f_dict['top'] = []
only_f_dict['top_vs'] = []
for x in approach_method_metrics:
    only_f_dict[x] = []
    only_f_dict[x + '_vs'] = []

for f_id in range(0, len(only_f_bops)):
    bop = only_f_bops[f_id]
    bop_id = f_top_dict['bop'].index(bop)
    is_in_vs = True if bop in m_all_dict['bop'] else False
    if is_in_vs:
        bop_id_vs = m_all_dict['bop'].index(bop)
    else:
        bop_id_vs = -1

    only_f_dict['bop'].append(bop)
    only_f_dict['gene'].append(f_top_dict['gene'][bop_id])
    only_f_dict['top'].append(bop_id)
    if is_in_vs:
        only_f_dict['top_vs'].append(bop_id_vs)
    else:
        only_f_dict['top_vs'].append('None')
    for x in approach_method_metrics:
        only_f_dict[x].append(f_top_dict[x][bop_id])
        if is_in_vs:
            only_f_dict[x + '_vs'].append(m_all_dict[x][bop_id_vs])
        else:
            only_f_dict[x + '_vs'].append('None')

only_f_order = np.argsort(only_f_dict['top'])
only_f_dict['bop'] = list(np.array(only_f_dict['bop'])[only_f_order])
only_f_dict['gene'] = list(np.array(only_f_dict['gene'])[only_f_order])
only_f_dict['top'] = list(np.array(only_f_dict['top'])[only_f_order])
only_f_dict['top_vs'] = list(np.array(only_f_dict['top_vs'])[only_f_order])

for x in approach_method_metrics:
    only_f_dict[x] = list(np.array(only_f_dict[x])[only_f_order])
    only_f_dict[x + '_vs'] = list(np.array(only_f_dict[x + '_vs'])[only_f_order])

only_f_df = pd.DataFrame(only_f_dict)

# Table for male-only
only_m_dict = {}
only_m_dict['bop'] = []
only_m_dict['gene'] = []
only_m_dict['top'] = []
only_m_dict['top_vs'] = []
for x in approach_method_metrics:
    only_m_dict[x] = []
    only_m_dict[x + '_vs'] = []

for f_id in range(0, len(only_m_bops)):
    bop = only_m_bops[f_id]
    bop_id = m_top_dict['bop'].index(bop)

    is_in_vs = True if bop in f_all_dict['bop'] else False
    if is_in_vs:
        bop_id_vs = f_all_dict['bop'].index(bop)
    else:
        bop_id_vs = -1

    only_m_dict['bop'].append(bop)
    only_m_dict['gene'].append(m_top_dict['gene'][bop_id])
    only_m_dict['top'].append(bop_id)
    if is_in_vs:
        only_m_dict['top_vs'].append(bop_id_vs)
    else:
        only_m_dict['top_vs'].append('None')
    for x in approach_method_metrics:
        only_m_dict[x].append(m_top_dict[x][bop_id])
        if is_in_vs:
            only_m_dict[x + '_vs'].append(f_all_dict[x][bop_id_vs])
        else:
            only_m_dict[x + '_vs'].append('None')

only_m_order = np.argsort(only_m_dict['top'])
only_m_dict['bop'] = list(np.array(only_m_dict['bop'])[only_m_order])
only_m_dict['gene'] = list(np.array(only_m_dict['gene'])[only_m_order])
only_m_dict['top'] = list(np.array(only_m_dict['top'])[only_m_order])
only_m_dict['top_vs'] = list(np.array(only_m_dict['top_vs'])[only_m_order])

for x in approach_method_metrics:
    only_m_dict[x] = list(np.array(only_m_dict[x])[only_m_order])
    only_m_dict[x + '_vs'] = list(np.array(only_m_dict[x + '_vs'])[only_m_order])

only_m_df = pd.DataFrame(only_m_dict)

# Table for intersection
i_dict = {}
i_dict['bop'] = []
i_dict['gene'] = []
i_dict['top_f'] = []
i_dict['top_m'] = []
for x in approach_method_metrics:
    i_dict[x + '_f'] = []
    i_dict[x + '_m'] = []

for f_id in range(0, len(i_bops)):
    bop = i_bops[f_id]
    bop_id_f = f_top_dict['bop'].index(bop)
    bop_id_m = m_top_dict['bop'].index(bop)

    i_dict['bop'].append(bop)
    i_dict['gene'].append(f_top_dict['gene'][bop_id_f])
    i_dict['top_f'].append(bop_id_f)
    i_dict['top_m'].append(bop_id_m)
    for x in approach_method_metrics:
        i_dict[x + '_f'].append(f_top_dict[x][bop_id_f])
        i_dict[x + '_m'].append(m_top_dict[x][bop_id_m])

i_order = np.argsort(i_dict['top_f'])
i_dict['bop'] = list(np.array(i_dict['bop'])[i_order])
i_dict['gene'] = list(np.array(i_dict['gene'])[i_order])
i_dict['top_f'] = list(np.array(i_dict['top_f'])[i_order])
i_dict['top_m'] = list(np.array(i_dict['top_m'])[i_order])

for x in approach_method_metrics:
    i_dict[x + '_f'] = list(np.array(i_dict[x + '_f'])[i_order])
    i_dict[x + '_m'] = list(np.array(i_dict[x + '_m'])[i_order])

i_df = pd.DataFrame(i_dict)

writer = pd.ExcelWriter(approach_method.value + '.xlsx', engine='xlsxwriter')
only_f_df.to_excel(writer, index=False, sheet_name='only_f')
only_m_df.to_excel(writer, index=False, sheet_name='only_m')
i_df.to_excel(writer, index=False, sheet_name='i')
writer.save()
