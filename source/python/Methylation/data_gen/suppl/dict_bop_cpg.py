from annotations.bop import get_dict_bop_cpgs
from config.config import *

def save_dict_bop_cpg(config):
    dict_bop_cpgs = get_dict_bop_cpgs(config)

    str_list = []
    for bop in dict_bop_cpgs:
        bop_cpgs = dict_bop_cpgs.get(bop)
        curr_str = bop
        for bop_cpg in bop_cpgs:
            curr_str += ' ' + bop_cpg
        str_list.append(curr_str)

    fn = get_bop_data_path(config, 'dict_bop_cpg.txt')
    np.savetxt(fn, str_list, fmt="%s")


db = DataBaseType.GSE52588
config = Config(
    db=db
)

classes = [ClassType.any, ClassType.class_a, ClassType.class_b, ClassType.class_c, ClassType.class_d]

for class_type in classes:
    print('class: ' + class_type.value)
    config.class_type = class_type
    save_dict_bop_cpg(config)
