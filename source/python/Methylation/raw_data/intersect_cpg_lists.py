from infrastructure.save.features import *
import pandas as pd

def intersect_xlsx(files_names):

    keys = []

    data = dict()
    for file_name in files_names:
        key = file_name[:-5]
        keys.append(key)
        data[key] = []
        df = pd.read_excel(file_name)
        data[key].append(list(df.values))


    data_intersection = []
    data_column_names = []
    data_column_names.append('cpg')
    for key in keys:
        data_column_names.append(key[:8])
    for target_line in data[keys[0]][0]:
        cpg = target_line[0]
        for key in keys[1:]:
            for item in data[key][0]:
                if cpg == item[0]:
                    item = list(item)
                    item.append(target_line[1])
                    data_intersection.append(item)

    fn = 'common' + keys[0][8:] + '.xlsx'
    df = pd.DataFrame(data_intersection)
    df.columns = data_column_names
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

files_names = ['GSE40279_by_slope_intersection.xlsx',
               'GSE87571_by_slope_intersection.xlsx']
intersect_xlsx(files_names)