from config.config import *
from config.types.annotations import GeneDataType, GeoType
from config.types.attributes.common import Gender
from infrastructure.load.top import *
import xlsxwriter


num_top = 500

db = DataBase.GSE40279
dt = DataType.gene
approach = Approach.top
scenario = Scenario.approach
geo = GeoType.islands_shores
genders = [Gender.any, Gender.M, Gender.F]

fn = 'krivonosov_1_genes.txt'
gene_top_origin = load_top_gene_names_by_article(Config(read_only=True), fn)

column = 0
places = {}
configs = {}

approach_gd = GeneDataType.mean
approach_methods = [Method.anova, Method.enet, Method.linreg, Method.spearman]
for gender in genders:
    print('gender: ' + gender.value)
    for method in approach_methods:
        print('\t' + 'method: ' + method.value)

        config = Config(
            read_only=True,
            db=db,
            dt=dt,
            approach=approach,
            scenario=scenario,
            approach_method=method,
            gender=gender,
            approach_gd=approach_gd,
            geo=geo
        )

        gene_top = load_top_gene_names(config, num_top)[0:num_top]

        curr_places = []
        for gene in gene_top_origin:
            if gene in gene_top:
                curr_places.append(gene_top.index(gene))
            else:
                curr_places.append('nan')

        config_name = 'data_type: ' + approach_gd.value + '; ' + \
                      'method: ' + method.value + '; ' + \
                      'gender: ' + gender.value

        places[column] = curr_places
        configs[column] = config_name

        column += 1

approach_gd = GeneDataType.from_cpg
approach_methods = [Method.anova, Method.enet, Method.linreg, Method.spearman]
for gender in genders:
    print('gender: ' + gender.value)
    for method in approach_methods:
        print('\t' + 'method: ' + method.value)

        config = Config(
            read_only=True,
            db=db,
            dt=dt,
            approach=approach,
            scenario=scenario,
            approach_method=method,
            gender=gender,
            approach_gd=approach_gd,
            geo=geo
        )

        gene_top = load_top_gene_names(config, num_top)[0:num_top]

        curr_places = []
        for gene in gene_top_origin:
            if gene in gene_top:
                curr_places.append(gene_top.index(gene))
            else:
                curr_places.append('nan')

        config_name = 'data_type: ' + approach_gd.value + '; ' + \
                      'method: ' + method.value + '; ' + \
                      'gender: ' + gender.value

        places[column] = curr_places
        configs[column] = config_name

        column += 1

approach_gd = GeneDataType.from_bop
approach_methods = [Method.manova]
for gender in genders:
    print('gender: ' + gender.value)
    for method in approach_methods:
        print('\t' + 'method: ' + method.value)

        config = Config(
            read_only=True,
            db=db,
            dt=dt,
            approach=approach,
            scenario=scenario,
            approach_method=method,
            gender=gender,
            approach_gd=approach_gd,
            geo=geo
        )

        gene_top = load_top_gene_names(config, num_top)[0:num_top]

        curr_places = []
        for gene in gene_top_origin:
            if gene in gene_top:
                curr_places.append(gene_top.index(gene))
            else:
                curr_places.append('nan')

        config_name = 'data_type: ' + approach_gd.value + '; ' + \
                      'method: ' + method.value + '; ' + \
                      'gender: ' + gender.value

        places[column] = curr_places
        configs[column] = config_name

        column += 1


places = list(places.values())
places = list(map(list, zip(*places)))

print(configs.keys())
first_row = list(configs.values())
table = [first_row]

for row in places:
    table.append(row)\

table = list(map(list, zip(*table)))

fn = 'match.xlsx'
fn = get_path(Config(db=db, read_only=True), fn)
workbook = xlsxwriter.Workbook(fn)
worksheet = workbook.add_worksheet()
row = 0
for col, data in enumerate(table):
    worksheet.write_column(row, col, data)
workbook.close()
