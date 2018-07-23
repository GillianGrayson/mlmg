from config.types import DataBaseType
from infrastructure.file_system import get_path
import numpy as np


def load_attributes(config):
    fn = 'attributes.txt'
    fn = get_path(config, fn)
    f = open(fn)

    key_line = f.readline()
    keys = key_line.split(' ')
    attributes = {}
    for key_id in range(0, len(keys)):
        keys[key_id] = keys[key_id].rstrip()
        attributes[keys[key_id]] = []

    for line in f:
        vals = line.split(' ')
        for key_id in range(0, len(keys)):
            attributes[keys[key_id]].append(vals[key_id].rstrip())

    return attributes

def get_main_attributes(config):
    indexes = config.indexes
    ma = []
    if config.db is DataBaseType.GSE40279:
        ma = list(map(int, list(np.array(config.attributes['age'])[indexes])))
    elif config.db is DataBaseType.GSE52588:
        ma = config.attributes['type']
    return ma