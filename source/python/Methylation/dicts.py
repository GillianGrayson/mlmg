from infrastructure.file_system import *
from geo import *

def get_dicts(config):
    fs_type = config.fs_type
    db_type = config.db_type
    geo_type = config.geo_type

    target_geo = []
    if geo_type is GeoType.islands:
        target_geo.append('Island')
    elif geo_type is GeoType.shores:
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')
    elif geo_type is GeoType.shores_s:
        target_geo.append('S_Shore')
    elif geo_type is GeoType.shores_n:
        target_geo.append('N_Shore')
    elif geo_type is GeoType.islands_shores:
        target_geo.append('Island')
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')

    fn = 'cpg.txt'
    cpg = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            cpg.append(line)

    fn = 'gene.txt'
    gene = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            gene.append(line)

    fn = 'chr.txt'
    chr = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            chr.append(line)

    fn = 'cpg_type.txt'
    cpg_type = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            cpg_type.append(line)

    fn = 'map_info.txt'
    map_info = []
    full_path = get_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            map_info.append(line)

    on_sex_chr = 0
    dict_cpg_gene = {}
    dict_cpg_map = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_cpg_type = cpg_type[i].rstrip()
        curr_map_info = map_info[i].rstrip()

        is_target_geo = False
        if len(target_geo) == 0:
            is_target_geo = True
        else:
            if curr_cpg_type in target_geo:
                is_target_geo = True

        if len(curr_cpg) > 2:
            if curr_cpg[0:2] == 'cg' or curr_cpg[0:2] == 'rs' or curr_cpg[0:2] == 'ch':
                if curr_chr != 'X' and curr_chr != 'Y':
                    if is_target_geo:
                        if len(curr_gene) > 0:
                            all_genes = list(set(curr_gene.split(';')))
                            dict_cpg_gene[curr_cpg] = all_genes
                            dict_cpg_map[curr_cpg] = curr_map_info
                else:
                    on_sex_chr += 1

    print('on_sex_chr: ' + str(on_sex_chr))

    return dict_cpg_gene, dict_cpg_map
