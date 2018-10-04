from infrastructure.path.path import *
from config.types import *
import os.path
import pickle


def load_annotations(config):
    fn_txt = 'annotations.txt'
    fn_pkl = 'annotations.pkl'
    fn_txt = get_path(config, fn_txt)
    fn_pkl = get_path(config, fn_pkl)

    is_pkl = os.path.isfile(fn_pkl)

    if is_pkl:
        f = open(fn_pkl, 'rb')
        annotations = pickle.load(f)
        f.close()
    else:
        f = open(fn_txt)

        key_line = f.readline()
        keys = key_line.split('\t')
        annotations = {}
        for key_id in range(0, len(keys)):
            keys[key_id] = keys[key_id].rstrip()
            annotations[keys[key_id]] = []

        for line in f:
            vals = line.split('\t')
            for key_id in range(0, len(keys)):
                key = keys[key_id]
                if key == Annotation.cpg.value:
                    annotations[key].append(vals[key_id].rstrip())
                elif key == Annotation.chr.value:
                    annotations[key].append(vals[key_id].rstrip())
                elif key == Annotation.map_info.value:
                    if vals[key_id].rstrip().isdigit():
                        annotations[key].append(int(vals[key_id].rstrip()))
                    else:
                        annotations[key].append(0)
                elif key == Annotation.gene.value:
                    annotations[key].append(vals[key_id].rstrip())
                elif key == Annotation.class_type.value:
                    annotations[key].append(vals[key_id].rstrip())
                elif key == Annotation.geo.value:
                    annotations[key].append(vals[key_id].rstrip())
                elif key == Annotation.bop.value:
                    annotations[key].append(vals[key_id].rstrip())

        f.close()

    if not is_pkl:
        f = open(fn_pkl, 'wb')
        pickle.dump(annotations, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    return annotations
