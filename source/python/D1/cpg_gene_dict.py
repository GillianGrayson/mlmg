from Infrastructure.file_system import *

def get_dict_cpg_gene(type):

    fn = 'cpg.txt'
    cpg = []
    full_path = get_full_path(type, fn)
    with open(full_path) as f:
        for line in f:
            cpg.append(line)

    fn = 'gene.txt'
    gene = []
    full_path = get_full_path(type, fn)
    with open(full_path) as f:
        for line in f:
            gene.append(line)

    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()

        if len(curr_cpg) > 2:
            if curr_cpg[0:2] == 'cg':
                if len(curr_gene) > 0:
                    all_genes = list(set(curr_gene.split(';')))
                    dict_cpg_gene[curr_cpg] = all_genes

    return dict_cpg_gene



