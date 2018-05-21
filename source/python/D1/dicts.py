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

    fn = 'chr.txt'
    chr = []
    full_path = get_full_path(type, fn)
    with open(full_path) as f:
        for line in f:
            chr.append(line)

    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()

        if len(curr_cpg) > 2:
            if curr_cpg[0:2] == 'cg':
                if curr_chr != 'X' and curr_chr != 'Y':
                    if len(curr_gene) > 0:
                        all_genes = list(set(curr_gene.split(';')))
                        dict_cpg_gene[curr_cpg] = all_genes

    return dict_cpg_gene

def get_dict_cpg_gene_shore(type):

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

    fn = 'chr.txt'
    chr = []
    full_path = get_full_path(type, fn)
    with open(full_path) as f:
        for line in f:
            chr.append(line)

    fn = 'cpg_type.txt'
    cpg_type = []
    full_path = get_full_path(type, fn)
    with open(full_path) as f:
        for line in f:
            cpg_type.append(line)

    dict_cpg_gene = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()
        curr_chr = chr[i].rstrip()
        curr_cpg_type = cpg_type[i].rstrip()

        if len(curr_cpg) > 2:
            if curr_cpg[0:2] == 'cg':
                if curr_cpg_type == 'N_Shore' or curr_cpg_type == 'S_Shore':
                    if curr_chr != 'X' and curr_chr != 'Y':
                        if len(curr_gene) > 0:
                            all_genes = list(set(curr_gene.split(';')))
                            dict_cpg_gene[curr_cpg] = all_genes

    return dict_cpg_gene

def get_dict_gene_cpg(type):

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

    dict_gene_cpg = {}
    for i in range(0, len(cpg)):

        curr_cpg = cpg[i].rstrip()
        curr_gene = gene[i].rstrip()

        if len(curr_cpg) > 2:
            if curr_cpg[0:2] == 'cg':
                if len(curr_gene) > 0:
                    all_genes = list(set(curr_gene.split(';')))
                    for key in all_genes:
                        if key in dict_gene_cpg:
                            dict_gene_cpg[key].append(curr_cpg)
                        else:
                            dict_gene_cpg[key] = [curr_cpg]

    for key in dict_gene_cpg:
        dict_gene_cpg[key] = list(set(dict_gene_cpg[key]))

    return dict_gene_cpg
