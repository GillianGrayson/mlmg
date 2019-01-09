from config.config import *
from annotations.bop import *
from config.types.attributes.common import Disease, Gender
from infrastructure.load.top import *
from config.types.experiments.method import *
import pandas as pd
from infrastructure.path.path import get_origin_path

def gene_from_bop(config, name):

    fn = name + '.txt'
    file = open(fn)
    bops = file.read().splitlines()
    file.close()
    dict_bop_genes = get_dict_bop_genes(config)

    genes = []

    for bop in bops:
        if bop in dict_bop_genes:
            curr_genes = dict_bop_genes[bop]
            for gene in curr_genes:
                if gene not in genes:
                    genes.append(gene)

    save_dict = {'name': genes}

    fn = name + '.xlsx'
    df = pd.DataFrame(save_dict)
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()


data_base = DataBase.GSE87571
data_type = DataType.gene

cross_reactive = CrossReactiveType.cross_reactive_included
snp = SNPType.snp_included

chromosome_type = ChromosomeType.non_gender

class_type = ClassType.class_ab
dna_region = DNARegionType.genic

geo_type = GeoType.islands_shores
gene_data_type = GeneDataType.mean

disease = Disease.any
gender = Gender.any

scenario = Scenario.approach
approach = Approach.top
method = Method.linreg_ols

is_clustering = False


config = Config(
    read_only=False,

    data_base=data_base,

    data_type=data_type,

    cross_reactive=cross_reactive,
    snp=snp,

    chromosome_type=chromosome_type,

    class_type=class_type,

    dna_region=dna_region,

    geo_type=geo_type,
    gene_data_type=gene_data_type,

    scenario=scenario,
    approach=approach,
    method=method,

    disease=disease,
    gender=gender,

    is_clustering=is_clustering
)

gene_from_bop(config, 'bop_F')

