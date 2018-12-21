from infrastructure.path.path import get_path
import os.path


def load_excluded(config):
    fn_txt = 'snp_cluster.txt'
    fn_txt = get_path(config, fn_txt)

    excluded = []

    if os.path.isfile(fn_txt):

        f = open(fn_txt)

        for line in f:
            excluded.append(line.rstrip())
        f.close()


    return excluded
