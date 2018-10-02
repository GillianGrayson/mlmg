import itertools
import os.path


def create_fs(categories_vals, categories_names, init_path):
    all = list(itertools.product(*categories_vals))
    for dir_list in all:
        path = init_path
        for id in range(0, len(dir_list)):
            create_fs_description_file(path + '/' + categories_names[id] + '.txt', categories_vals[id])
            path += '/' + dir_list[id]
            if not os.path.exists(path):
                os.makedirs(path)


def create_fs_description_file(fn, instances):
    f = open(fn, 'w')
    for x in instances:
        f.write(x + '\n')
    f.close()