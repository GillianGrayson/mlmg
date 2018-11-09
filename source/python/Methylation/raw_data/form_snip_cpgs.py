from config.config import *
from infrastructure.save.features import *

def form_snip_cpgs(config):
    fn_txt = 'annotations.txt'
    fn_txt = get_path(config, fn_txt)

    target_key = 'Probe_SNPs_10'

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
            elif key == target_key:
                annotations[key].append(vals[key_id].rstrip())
    f.close()

    ids = [id for id in range(0, len(annotations[target_key])) if annotations[target_key][id] != '']

    SNPs_cpgs = list(np.array(annotations[Annotation.cpg.value])[ids])

    fn = get_path(config,  target_key + '.txt')
    save_features(fn, [SNPs_cpgs])


config = Config(data_base=DataBase.GSE87571)
form_snip_cpgs(config)