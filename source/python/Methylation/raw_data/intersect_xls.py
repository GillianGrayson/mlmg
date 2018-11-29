import pandas as pd
import os

def intersect_xls(files_names, files_pathes):

    keys = []

    data = dict()
    for file_id in range(0, len(files_names)):
        file_name = files_names[file_id]
        file_path = files_pathes[file_id]
        key = file_name[:-5]
        keys.append(key)
        data[key] = []
        df = pd.read_excel(file_path)
        data[key] = list(df.names)

    data_intersection = []
    data_column_names = []
    data_column_names.append('cpg')
    for target_line in data[keys[0]]:
        cpg = target_line
        cpg_check = [True]
        for key in keys[1:]:
            is_cpg_found = False
            for item in data[key]:
                if cpg == item:
                    is_cpg_found = True
            cpg_check.append(is_cpg_found)
        if all(cpg_check):
            data_intersection.append(cpg)

    fn = 'common.xlsx'
    df = pd.DataFrame(data_intersection)
    df.columns = data_column_names
    writer = pd.ExcelWriter(fn, engine='xlsxwriter')
    df.to_excel(writer, index=False)
    writer.save()

path = os.path.dirname(os.path.abspath(__file__))
dir_name = 'xls'
path = path + '\\' + dir_name

files_names = os.listdir(path)
files_pathes = [path + '\\' + file_name for file_name in files_names]
intersect_xls(files_names, files_pathes)