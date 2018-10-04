from config.config import *
from infrastructure.load.top import *
from config.method import *
import pandas as pd

num_top = 1000

data_base_type = DataBase.GSE40279

scenario = Scenario.approach
approach = Approach.top
disease = Disease.any

bop_data_type = DataType.bop
bop_approach_method = Method.manova
bop_metrics = get_method_metrics(bop_approach_method)

gene_data_type = DataType.gene
gene_approach_method = Method.linreg
gene_metrics = get_method_metrics(gene_approach_method)
approach_gene_data_type = GeneDataType.mean
geo_type = GeoType.islands_shores

bop_config_f = Config(
    read_only=True,
    db=data_base_type,
    dt=bop_data_type,
    scenario=scenario,
    approach=approach,
    approach_method=bop_approach_method,
    gender=Gender.F,
    disease=disease
)

bop_config_m = Config(
    read_only=True,
    db=data_base_type,
    dt=bop_data_type,
    scenario=scenario,
    approach=approach,
    approach_method=bop_approach_method,
    gender=Gender.M,
    disease=disease
)

config_f = Config(
    read_only=True,
    db=data_base_type,
    dt=gene_data_type,
    scenario=scenario,
    approach=approach,
    approach_method=gene_approach_method,
    gender=Gender.F,
    disease=disease,
    approach_gd=approach_gene_data_type,
    geo=geo_type
)

config_m = Config(
    read_only=True,
    db=data_base_type,
    dt=gene_data_type,
    scenario=scenario,
    approach=approach,
    approach_method=gene_approach_method,
    gender=Gender.M,
    disease=disease,
    approach_gd=approach_gene_data_type,
    geo=geo_type
)

bop_keys = ['bop', 'gene'] + bop_metrics
gene_keys = ['gene'] + gene_metrics

bop_f_top_dict = load_top_dict(bop_config_f, bop_keys, num_top)
bop_f_all_dict = load_top_dict(bop_config_f, bop_keys)
bop_m_top_dict = load_top_dict(bop_config_m, bop_keys, num_top)
bop_m_all_dict = load_top_dict(bop_config_m, bop_keys)
gene_f_top_dict = load_top_dict(config_f, gene_keys, num_top)
gene_f_all_dict = load_top_dict(config_f, gene_keys)
gene_m_top_dict = load_top_dict(config_m, gene_keys, num_top)
gene_m_all_dict = load_top_dict(config_m, gene_keys)

bop_f_genes_raw = bop_f_top_dict['gene']
bop_f_genes = []
for genes in bop_f_genes_raw:
    gene_list = genes.split(';')
    for gene in gene_list:
        if gene not in bop_f_genes:
            bop_f_genes.append(gene)
bop_m_genes_raw = bop_m_top_dict['gene']
bop_m_genes = []
for genes in bop_m_genes_raw:
    gene_list = genes.split(';')
    for gene in gene_list:
        if gene not in bop_m_genes:
            bop_m_genes.append(gene)
bop_i_genes = list(set(bop_f_genes).intersection(bop_m_genes))
bop_only_f_genes = list(set(bop_f_genes) - set(bop_i_genes))
bop_only_m_genes = list(set(bop_m_genes) - set(bop_i_genes))

gene_f_genes = gene_f_top_dict['gene']
gene_m_genes = gene_m_top_dict['gene']

gene_i_genes = list(set(gene_f_genes).intersection(gene_m_genes))
gene_only_f_genes = list(set(gene_f_genes) - set(gene_i_genes))
gene_only_m_genes = list(set(gene_m_genes) - set(gene_i_genes))

gene_vs_bop_i_genes = list(set(bop_i_genes).intersection(gene_i_genes))

a = 1



