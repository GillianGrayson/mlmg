from infrastructure.save.features import save_features
from annotations.gene import get_dict_gene_cpg
from infrastructure.load.cpg_data import load_dict_cpg_data
from config.config import *


def save_inside_gene_moment(config):
    dict_gene_cpg = get_dict_gene_cpg(config)
    dict_cpg_data = load_dict_cpg_data(config)

    genes_opt = []
    cpgs_opt = []
    means_opt = []
    stds_opt = []
    mins_opt = []
    maxs_opt = []
    cpg_id = 0
    for gene, cpgs in dict_gene_cpg.items():
        for cpg in cpgs:
            if cpg in dict_cpg_data:
                cpg_data = dict_cpg_data.get(cpg)

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
        genes_opt,
        cpgs_opt,
        means_opt,
        stds_opt,
        mins_opt,
        maxs_opt
    ]
    fn = 'inside_gene.txt'
    fn = get_result_path(config, fn)
    save_features(fn, features)


data_base = DataBase.GSE87571
data_type = DataType.gene

cross_reactives = [CrossReactiveType.cross_reactive_excluded_weak]
snps = [SNPType.snp_excluded_weak]

chromosome_type = ChromosomeType.non_gender

geo_types = [GeoType.any]
gene_data_type = GeneDataType.mean

scenario = Scenario.approach
approach = Approach.inside_gene
method = Method.moment

disease = Disease.any
genders = [Gender.F, Gender.M, Gender.any]

for cross_reactive in cross_reactives:
    print(cross_reactive.value)
    for snp in snps:
        print(snp.value)
        for gender in genders:
            print('\t' + gender.value)
            for geo_type in geo_types:
                print('\t\t' + geo_type.value)

                config = Config(
                    data_base=data_base,
                    data_type=data_type,

                    cross_reactive=cross_reactive,
                    snp=snp,

                    chromosome_type=chromosome_type,

                    geo_type=geo_type,
                    gene_data_type=gene_data_type,

                    scenario=scenario,
                    approach=approach,
                    method=method,

                    disease=disease,
                    gender=gender
                )

                save_inside_gene_moment(config)
