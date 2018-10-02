from config.types import *
from infrastructure.path.path import get_path
import numpy as np
import os.path
import pickle


def load_attributes(config):
    fn_txt = 'attributes.txt'
    fn_pkl = 'attributes.pkl'
    fn_txt = get_path(config, fn_txt)
    fn_pkl = get_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)
    if is_pkl:
        f = open(fn_pkl, 'rb')
        attributes = pickle.load(f)
        f.close()
    else:

        f = open(fn_txt)

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

        f = open(fn_pkl, 'wb')
        pickle.dump(attributes, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    return attributes

def get_attributes(config, attribute=Attribute.age):
    indexes = config.indexes
    atr = []
    if attribute is Attribute.age:
        atr = list(map(int, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.gender:
        atr = list(map(int, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.disease:
        atr = list(map(str, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.group:
        atr = list(map(int, list(np.array(config.attributes[attribute.value])[indexes])))
    elif attribute is Attribute.batch:
        atr = list(map(int, list(np.array(config.attributes[attribute.value])[indexes])))

    return atr
