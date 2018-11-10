import pandas as pd

def intersect_xlsx(files_names):

    keys = []
    column_names = dict()
    data = dict()
    for file_name in files_names:
        key = file_name[:-5]
        keys.append(key)
        data[key] = []
        column_names[key] = []
        df = pd.read_excel(file_name)
        data[key].append(list(df.values))
        column_names[key].append(list(df.columns))


    data_intersection = []
    data_column_names = []
    for column_name in column_names[keys[1]][0]:
        data_column_names.append(column_name)
    for target_line in data[keys[0]][0]:
        cpg = target_line[0]
        for key in keys[1:]:
            for item in data[key][0]:
                if cpg == item[0]:
                    item = list(item)
                    data_intersection.append(item)

    fn = keys[1] + '_' + keys[0] + '.xlsx'
    df = pd.DataFrame(data_intersection)
    df.columns = data_column_names
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

files_names = ['by_condition.xlsx',
               'method(linreg_ols)_wo_cross_reactive.xlsx']
intersect_xlsx(files_names)