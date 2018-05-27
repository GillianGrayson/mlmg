from enum import Enum

class FSType(Enum):
    local_big = 'E:/Work/mlmg/data'
    local_msi = 'D:/Work/mlmg/data'
    unn = '/common/home/yusipov_i/Work/mlmg/data'
    mpipks = '/data/biophys/yusipov/mlmg/data'

class DataBaseType(Enum):
    GSE40279 = 'GSE40279'
    GSE52588 = 'GSE52588'

def get_root(fs_type):
    root = fs_type.value
    return root

def get_full_path(fs_type, db_type, file_name):
    root = fs_type.value
    db = db_type.value
    path = root + '/' + db + '/' + file_name
    return path