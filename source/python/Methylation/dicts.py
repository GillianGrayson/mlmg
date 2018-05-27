from Infrastructure.file_system import *
from gen_files.geo import *

def get_dict_cpg_gene(fs_type, db_type, geo_type):

    target_geo = []
    if geo_type is GeoType.islands:
        target_geo.append('Island')
    elif geo_type is GeoType.shores:
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')
    elif geo_type is GeoType.islands_shores:
        target_geo.append('Island')
        target_geo.append('N_Shore')
        target_geo.append('S_Shore')

    fn = 'cpg.txt'
    cpg = []
    full_path = get_full_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            cpg.append(line)

    fn = 'gene.txt'
    gene = []
    full_path = get_full_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            gene.append(line)

    fn = 'chr.txt'
    chr = []
    full_path = get_full_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            chr.append(line)

    fn = 'cpg_type.txt'
    cpg_type = []
    full_path = get_full_path(fs_type, db_type, fn)
    with open(full_path) as f:
        for line in f:
            cpg_type.append(line)

    on_sex_chr = 0
    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_cpg_type = cpg_type[i].rstrip()

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
                else:
                    on_sex_chr += 1

    print('on_sex_chr: ' + str(on_sex_chr))

    return dict_cpg_gene
