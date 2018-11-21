from enum import Enum


class FSType(Enum):
    local_big = 'E:/YandexDisk/Work/mlmg/data'
    local_msi = 'D:/YandexDisk/Work/mlmg/data'
    local_ab = 'C:/Users/User/YandexDisk/mlmg/data'
    unn = '/common/home/yusipov_i/Work/mlmg/data'
    mpipks = '/data/biophys/yusipov/mlmg/data'
    unn_ab = '/common/home/kalyakulina_a/Work/mlmg/data'


class DataBase(Enum):
    GSE40279 = 'GSE40279'
    GSE52588 = 'GSE52588'
    GSE30870 = 'GSE30870'
    GSE61256 = 'GSE61256'
    GSE63347 = 'GSE63347'
    GSE52588_TEST = 'GSE52588_TEST'
    GSE87571 = 'GSE87571'
    data_base_versus = 'versus'


class DataType(Enum):
    cpg = 'cpg'
    gene = 'gene'
    bop = 'bop'
    versus = 'versus'


class InfoType(Enum):
    result = 'result'
    param = 'param'
    data = 'data'


