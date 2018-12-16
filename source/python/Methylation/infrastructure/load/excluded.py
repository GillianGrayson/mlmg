from infrastructure.path.path import get_path


def load_excluded(config):
    fn_txt = 'snp_cluster.txt'
    fn_txt = get_path(config, fn_txt)

    f = open(fn_txt)

    excluded = []
    for line in f:
        excluded.append(line.rstrip())
    f.close()

    return excluded
