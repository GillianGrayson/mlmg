from infrastructure.path import *
from config.types import *


def load_annotations(config):
    fn = 'annotations.txt'
    fn = get_path(config, fn)
    f = open(fn)

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

    return annotations
