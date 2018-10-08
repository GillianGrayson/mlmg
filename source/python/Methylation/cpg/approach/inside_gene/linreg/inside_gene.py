from infrastructure.save.features import save_features
from annotations.gene import get_dict_gene_cpg
from infrastructure.load.cpg_data import load_dict_cpg_data
from config.config import *
from scipy import stats


def save_inside_gene_linreg(config):
    attributes = get_attributes(config)
    dict_gene_cpg = get_dict_gene_cpg(config)
    dict_cpg_data = load_dict_cpg_data(config)

    genes_opt = []
    cpgs_opt = []
    r_values_opt = []
    p_values_opt = []
    slopes_opt = []
    intercepts_opt = []
    cpg_id = 0
    for gene, cpgs in dict_gene_cpg.items():
        for cpg in cpgs:
            cpg_data = dict_cpg_data.get(cpg)
            slope, intercept, r_value, p_value, std_err = stats.linregress(attributes, cpg_data)
            genes_opt.append(gene)
            cpgs_opt.append(cpg)
            r_values_opt.append(r_value)
            p_values_opt.append(p_value)
            slopes_opt.append(slope)
            intercepts_opt.append(intercept)
            if cpg_id % config.print_rate == 0:
                print('cpg_id: ' + str(cpg_id))
            cpg_id += 1

    features = [
        genes_opt,
        cpgs_opt,
        r_values_opt,
        p_values_opt,
        slopes_opt,
        intercepts_opt
    ]
    fn = 'inside_gene.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)


data_base = DataBase.GSE87571
data_type = DataType.cpg

chromosome_type = ChromosomeTypes.non_gender

dna_region = DNARegion.genic

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

scenario = Scenario.approach
approach = Approach.inside_gene
method = Method.linreg

for gender in genders:
    print('gender: ' + gender.value)

    config = Config(
        data_base=data_base,
        data_type=data_type,

        chromosome_type=chromosome_type,

        dna_region=dna_region,

        disease=disease,
        gender=gender,

        scenario=scenario,
        approach=approach,
        method=method
    )

    save_inside_gene_linreg(config)
