from config.types import *


def line_proc(config, line):
    line_list = []
    if config.db is DataBaseType.GSE40279:
        line_list = line.split('\t')
    elif config.db is DataBaseType.GSE52588:
        line_list = line.split('\t')
        for val_id in range(0, len(line_list)):
            line_list[val_id] = line_list[val_id].replace('"', '').rstrip()
    return line_list