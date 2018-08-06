import numpy as np
from infrastructure.path import *

def save_features(fn, features):
    features = np.array(features)
    shape = features.shape
    if len(shape) > 1:
        num_features = len(features)
        length = len(features[0])
        str_list = []
        for str_id in range(0, length):
            curr_str = str(features[0, str_id])
            for f_id in range(1, num_features):
                curr_str += ' ' + str(features[f_id, str_id])
            str_list.append(curr_str)
    else:
        length = len(features)
        str_list = []
        for str_id in range(0, length):
            curr_str = str(features[str_id])
            str_list.append(curr_str)
    np.savetxt(fn, str_list, fmt="%s")