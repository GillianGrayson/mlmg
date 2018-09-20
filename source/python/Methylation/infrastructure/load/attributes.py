from config.types import *
from infrastructure.path import get_path
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

    f.close()

    return attributes

def get_attributes(config, attribute=Attribute.age):
    indexes = config.indexes
    atr = []
    if attribute is Attribute.age:
        atr = list(map(int, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.gender:
        atr = list(map(str, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.disease:
        atr = list(map(str, list(np.array(config.attributes[attribute.value])[indexes])))

    return atr
