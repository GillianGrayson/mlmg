from infrastructure.save.features import save_features
from annotations.bop import get_dict_bop_cpgs, get_dict_bop_genes
from infrastructure.load.cpg_data import load_dict_cpg_data
from config.config import *


def save_inside_bop_moment(config):
    dict_bop_cpg = get_dict_bop_cpgs(config)
    dict_bop_gene = get_dict_bop_genes(config)
    dict_cpg_data = load_dict_cpg_data(config)

    bops_opt = []
    genes_opt = []
    cpgs_opt = []
    means_opt = []
    stds_opt = []
    mins_opt = []
    maxs_opt = []
    cpg_id = 0
    for bop, cpgs in dict_bop_cpg.items():
        if bop in dict_bop_gene:
            gene = ';'.join(dict_bop_gene.get(bop))
        else:
            gene = 'no_gene'
        for cpg in cpgs:
            if cpg in dict_cpg_data:
                cpg_data = dict_cpg_data.get(cpg)

                bops_opt.append(bop)
                genes_opt.append(gene)
                cpgs_opt.append(cpg)
                means_opt.append(np.mean(cpg_data))
                stds_opt.append(np.std(cpg_data))
                mins_opt.append(np.min(cpg_data))
                maxs_opt.append(np.max(cpg_data))

                if cpg_id % config.print_rate == 0:
                    print('cpg_id: ' + str(cpg_id))
                cpg_id += 1

    features = [
        bops_opt,
        genes_opt,
        cpgs_opt,
        means_opt,
        stds_opt,
        mins_opt,
        maxs_opt
    ]
    fn = 'inside_bop.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)


data_base = DataBase.GSE87571
data_type = DataType.bop

chromosome_type = ChromosomeType.non_gender

class_type = ClassType.class_ab

scenario = Scenario.approach
approach = Approach.inside_bop
method = Method.moment

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

for gender in genders:
    print('gender: ' + gender.value)

    config = Config(
        data_base=data_base,
        data_type=data_type,

        chromosome_type=chromosome_type,

        class_type=class_type,

        disease=disease,
        gender=gender,

        scenario=scenario,
        approach=approach,
        method=method
    )

    save_inside_bop_moment(config)
