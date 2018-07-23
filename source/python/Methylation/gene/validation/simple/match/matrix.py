from config.config import *
from infrastructure.load.top import *
from infrastructure.save.features import save_features
import xlsxwriter


db = DataBaseType.GSE40279
dt = DataType.gene
approach = Approach.top
scenario = Scenario.approach
approach_methods = [Method.anova, Method.enet, Method.linreg, Method.spearman]
approach_gd = GeneDataType.mean
geo = GeoType.islands_shores

config = Config(
    db=db,
    dt=dt,
    approach=approach,
    scenario=scenario,
    approach_method=Method.linreg,
    approach_gd=approach_gd,
    geo=geo
)

num_top = 100
fn = 'claudio2015.txt'
gene_top = load_top_gene_names_by_article(config, fn)
gene_top_dict = {}
gene_top_dict['claudio2015'] = gene_top

for method in approach_methods:
    config = Config(
        db=db,
        dt=dt,
        approach=approach,
        scenario=scenario,
        approach_method=method,
        approach_gd=approach_gd,
        geo=geo
    )
    gene_top = load_top_gene_names(config, num_top)[0:num_top]
    gene_top_dict[method.value] = gene_top

firsts = []
seconds = []
match_count = []
match_genes = []
first_row = list(gene_top_dict.keys())
first_row.insert(0, 'source')
table = [first_row]
for first in gene_top_dict:
    row = [first]
    for second in gene_top_dict:
        firsts.append(first)
        seconds.append(second)
        genes_first = gene_top_dict.get(first)
        genes_second = gene_top_dict.get(second)
        genes_common = []
        for gene in genes_first:
            if gene in genes_second:
                genes_common.append(gene)
        match_count.append(len(genes_common))
        match_genes.append(';'.join(genes_common))
        row.append(str(len(genes_common)))

    table.append(row)

config.scenario = Scenario.validation
config.validation = Validation.simple
config.validation_method = Method.match
config.validation_gd = config.approach_gd
fn = 'match_matrix.txt'
fn = get_result_path(config, fn)
save_features(fn, [firsts, seconds, match_count, match_genes])

fn = 'match_matrix.xlsx'
fn = get_result_path(config, fn)
workbook = xlsxwriter.Workbook(fn)
worksheet = workbook.add_worksheet()
row = 0
for col, data in enumerate(table):
    worksheet.write_column(row, col, data)
workbook.close()
