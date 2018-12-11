from config.config import *


def line_proc(config, line):
    line_list = []
    if config.data_base is DataBase.GSE40279:
        line_list = line.split('\t')
    elif config.data_base is DataBase.GSE52588 or config.data_base is DataBase.GSE52588_TEST:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    elif config.data_base is DataBase.GSE30870:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            if(line_list[val_id].rstrip() == ''):
                line_list[val_id] = config.miss_tag
    elif config.data_base is DataBase.GSE63347:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    elif config.data_base is DataBase.GSE87571:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    elif config.data_base is DataBase.liver:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    else:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    return line_list