from config.types import *


def line_proc(config, line):
    line_list = []
    if config.db is DataBaseType.GSE40279:
        line_list = line.split('\t')
    elif config.db is DataBaseType.GSE52588 or config.db is DataBaseType.GSE52588_TEST:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    elif config.db is DataBaseType.GSE30870:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            if(line_list[val_id].rstrip() == ''):
                line_list[val_id] = config.miss_tag
    elif config.db is DataBaseType.GSE63347:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    return line_list