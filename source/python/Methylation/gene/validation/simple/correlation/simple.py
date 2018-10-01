from config.config import *
from infrastructure.load.top import *
from annotations.gene import get_dict_gene_cpg
from infrastructure.load.cpg_data import load_cpg_data
from scipy import stats


num_top = 500

db = DataBaseType.GSE40279
dt = DataType.gene
approach = Approach.top
scenario = Scenario.approach
geo = GeoType.islands_shores
gender = Gender.F
approach_gd = GeneDataType.mean
approach_method = Method.linreg

num_genes = 10

config = Config(
    db=db,
    dt=dt,
    approach=approach,
    scenario=scenario,
    approach_method=approach_method,
    gender=gender,
    approach_gd=approach_gd,
    geo=geo
)

gene_top = load_top_gene_names(config, num_top)[0:num_top]

attributes = get_attributes(config)
dict_gene_cpg = get_dict_gene_cpg(config)
cpgs_names, cpg_vals = load_cpg_data(config)

for gene_id in range(0, num_genes):

    target_gene = gene_top[gene_id]
    target_cpgs = dict_gene_cpg.get(target_gene)

    curr_vals = 0
    prev_vals = 0

    print(target_gene)
    for cpg_id in range(0, len(target_cpgs)):
        prev_vals = curr_vals
        cpg = target_cpgs[cpg_id]
        curr_vals = cpg_vals[cpgs_names.index(cpg)]
        if cpg_id > 0:
            slope, intercept, r_value, p_value, std_err = stats.linregress(curr_vals, prev_vals)
            print('vs ' + str(r_value))
        slope, intercept, r_value, p_value, std_err = stats.linregress(curr_vals, attributes)
        print(cpg + ' ' + str(r_value))
    print('\n\n')
