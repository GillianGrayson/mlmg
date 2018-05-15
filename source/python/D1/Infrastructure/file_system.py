from enum import Enum

class FSType(Enum):
    local_big = 10
    local_msi = 20
    unn = 1
    mpipks = 2

def get_root(type):

    root = ''
    if type is FSType.unn:
        root = '/common/home/yusipov_i/Work/mlmg/data/D1'
    elif type is FSType.mpipks:
        root = '/data/biophys/yusipov/mlmg/data/D1'
    elif type is FSType.local_big:
        root = 'E:/Work/mlmg/data/D1'
    elif type is FSType.local_msi:
        root = 'D:/Work/mlmg/data/D1'

    return root

def get_full_path(type, file_name):
    root = get_root(type)
    path = root + '/' + file_name
    return path