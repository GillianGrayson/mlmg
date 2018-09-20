from config.config import *
from infrastructure.load.top import *

num_top = 500

method = Method.linreg

config_f = Config(
    read_only=True,
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
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
    db=DataBaseType.GSE40279,
    dt=DataType.gene,
    scenario=Scenario.approach,
    approach=Approach.top,
    approach_method=method,
    gender=Gender.M,
    disease=Disease.any,
    approach_gd=GeneDataType.mean,
    geo=GeoType.islands_shores
)

f_top_dict = load_top_gene_dict(config_f, num_top)
m_top_dict = load_top_gene_dict(config_m, num_top)

f_names = list(f_top_dict.keys())
m_names = list(m_top_dict.keys())

i_names = list(set(f_names).intersection(m_names))
only_f_names = list(set(f_names) - set(i_names))
only_m_names = list(set(m_names) - set(i_names))

only_f_ids = []
only_f_metrics = []
only_f_clusters = []
for f_id in range(0, len(only_f_names)):
    name = only_f_names[f_id]
    only_f_ids.append(f_top_dict[name][0])
    only_f_metrics.append(f_top_dict[name][1])
    only_f_clusters.append(f_top_dict[name][2])
only_f_order = np.argsort(only_f_ids)
only_f_ids_sorted = list(np.array(only_f_ids)[only_f_order])
only_f_names_sorted = list(np.array(only_f_names)[only_f_order])
only_f_metrics_sorted = list(np.array(only_f_metrics)[only_f_order])
only_f_clusters_sorted = list(np.array(only_f_clusters)[only_f_order])
strs = []
for f_id in range(0, len(only_f_names)):
    name = only_f_names_sorted[f_id]
    curr_str = name
    curr_str += (' ' + str(format(only_f_ids_sorted[f_id], 'd')))
    curr_str += (' ' + str(format(only_f_clusters_sorted[f_id], 'd')))
    curr_str += (' ' + str(format(only_f_metrics_sorted[f_id], '0.4f')))
    strs.append(curr_str)
fn = 'only_f.txt'
np.savetxt(fn, strs, fmt="%s")

only_m_ids = []
only_m_metrics = []
only_m_clusters = []
for m_id in range(0, len(only_m_names)):
    name = only_m_names[m_id]
    only_m_ids.append(m_top_dict[name][0])
    only_m_metrics.append(m_top_dict[name][1])
    only_m_clusters.append(m_top_dict[name][2])
only_m_order = np.argsort(only_m_ids)
only_m_ids_sorted = list(np.array(only_m_ids)[only_m_order])
only_m_names_sorted = list(np.array(only_m_names)[only_m_order])
only_m_metrics_sorted = list(np.array(only_m_metrics)[only_m_order])
only_m_clusters_sorted = list(np.array(only_m_clusters)[only_m_order])
strs = []
for m_id in range(0, len(only_m_names)):
    name = only_m_names_sorted[m_id]
    curr_str = name
    curr_str += (' ' + str(format(only_m_ids_sorted[m_id], 'd')))
    curr_str += (' ' + str(format(only_m_clusters_sorted[m_id], 'd')))
    curr_str += (' ' + str(format(only_m_metrics_sorted[m_id], '0.4f')))
    strs.append(curr_str)
fn = 'only_m.txt'
np.savetxt(fn, strs, fmt="%s")

i_f_ids = []
i_m_ids = []
i_f_clusters = []
i_m_clusters = []
i_diff_metrics = []
for i_id in range(0, len(i_names)):
    name = i_names[i_id]
    i_f_ids.append(f_top_dict[name][0])
    i_m_ids.append(m_top_dict[name][0])
    i_f_clusters.append(f_top_dict[name][2])
    i_m_clusters.append(m_top_dict[name][2])
    i_diff_metrics.append(f_top_dict[name][1] - m_top_dict[name][1])
i_order = np.argsort(i_f_ids)
i_names_sorted = list(np.array(i_names)[i_order])
i_f_ids_sorted = list(np.array(i_f_ids)[i_order])
i_m_ids_sorted = list(np.array(i_m_ids)[i_order])
i_f_clusters_sorted = list(np.array(i_f_clusters)[i_order])
i_m_clusters_sorted = list(np.array(i_m_clusters)[i_order])
i_diff_metrics_sorted = list(np.array(i_diff_metrics)[i_order])

strs = []
for i_id in range(0, len(i_names)):
    name = i_names_sorted[i_id]
    curr_str = name
    curr_str += (' ' + str(format(i_f_ids_sorted[i_id], 'd')))
    curr_str += (' ' + str(format(i_m_ids_sorted[i_id], 'd')))
    curr_str += (' ' + str(format(i_f_clusters_sorted[i_id], 'd')))
    curr_str += (' ' + str(format(i_m_clusters_sorted[i_id], 'd')))
    curr_str += (' ' + str(format(i_diff_metrics_sorted[i_id], '0.4f')))
    strs.append(curr_str)
fn = 'i.txt'
np.savetxt(fn, strs, fmt="%s")
