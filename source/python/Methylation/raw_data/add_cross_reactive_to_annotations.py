from config.config import *
from infrastructure.save.features import *

def add_cross_reactive_cpgs(config):
    cross_txt = 'cross_reactive_cpgs.txt'
    cross_txt = get_origin_path(config, cross_txt)
    cross_f = open(cross_txt)

    cross_reactive_cpgs = []
    for cpg in cross_f:
        cpg = cpg.replace('\n', '')
        cross_reactive_cpgs.append(cpg)

    cross_f.close()

    ann_txt = 'annotations.txt'
    ann_txt = get_path(config, ann_txt)
    ann_f = open(ann_txt)

    key_line = ann_f.readline()
    key_line = key_line.replace('\n', '')
    keys = key_line.split('\t')
    annotations = []
    lines = []

    for line in ann_f:
        line = line.replace('\n', '')
        lines.append(line)
        vals = line.split('\t')
        for key_id in range(0, len(keys)):
            key = keys[key_id]
            if key == Annotation.cpg.value:
                annotations.append(vals[key_id].rstrip())
    ann_f.close()

    key_line += '\tCROSS_R'

    for cpg_id in range(0, len(annotations)):
        print(cpg_id)
        cpg = annotations[cpg_id]
        if cpg in cross_reactive_cpgs:
            lines[cpg_id] += '\t1'
        else:
            lines[cpg_id] += '\t0'

    lines.insert(0, key_line)
    fn = get_path(config,  'annotations_mod.txt')
    save_features(fn, [lines])


config = Config(data_base=DataBase.GSE87571)
add_cross_reactive_cpgs(config)