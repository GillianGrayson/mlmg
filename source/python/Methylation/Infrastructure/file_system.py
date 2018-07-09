from enum import Enum
import socket

class FSType(Enum):
    local_big = 'E:/Work/mlmg/data'
    local_msi = 'D:/Work/mlmg/data'
    unn = '/common/home/yusipov_i/Work/mlmg/data'
    mpipks = '/data/biophys/yusipov/mlmg/data'

class DataBaseType(Enum):
    GSE40279 = 'GSE40279'
    GSE52588 = 'GSE52588'

class GeneDataType(Enum):
    mean = 'mean'
    mean_der = 'mean_der'
    mean_der_normed = 'mean_der_normed'
    from_cpg = 'cpg'

def get_root(fs_type):
    root = fs_type.value
    return root

def get_path(fs_type, db_type, file_name):
    root = fs_type.value
    db = db_type.value
    path = root + '/' + db + '/' + file_name
    return path

def get_gene_data_path(fs_type, db_type, gd_type, file_name):
    root = fs_type.value
    db = db_type.value
    gd = gd_type.value
    path = root + '/' + db + '/' + 'gene_data' + '/' + gd + '/' + file_name
    return path

def get_param_path(fs_type, db_type, file_name):
    root = fs_type.value
    db = db_type.value
    path = root + '/' + db + '/' + 'param' + '/' + file_name
    return path

def get_result_path(fs_type, db_type, file_name):
    root = fs_type.value
    db = db_type.value
    path = root + '/' + db + '/' + 'result' + '/' + file_name
    return path
