from config.types.attributes.cell_pop import CellPop
from infrastructure.path.path import get_path
import numpy as np


def load_cell_pop(config):
    fn = 'cell_pop.txt'
    fn = get_path(config, fn)
    f = open(fn)

    key_line = f.readline()
    keys = key_line.split('\t')
    cell_pop = {}
    for key_id in range(0, len(keys)):
        keys[key_id] = keys[key_id].rstrip()
        cell_pop[keys[key_id]] = []

    for line in f:
        vals = line.split('\t')
        for key_id in range(0, len(keys)):
            key = keys[key_id]
            if key == CellPop.sample_id.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.plasma_blast.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.cd8_p.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.cd8_naive.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.cd4_naive.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.cd8_t.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.cd4_t.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.nk.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.b_cell.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.mono.value:
                cell_pop[key].append(vals[key_id].rstrip())
            elif key == CellPop.gran.value:
                cell_pop[key].append(vals[key_id].rstrip())

    f.close()
    return cell_pop


def get_cell_pop(config, cell_pop):
    indexes = config.indexes
    pop = []
    for cell_pop_type in cell_pop:
        if cell_pop_type is CellPop.sample_id:
            pop = list(map(str, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.plasma_blast:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.cd8_p:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.cd8_naive:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.cd4_naive:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.cd8_t:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.cd4_t:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.nk:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.b_cell:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.mono:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))
        elif cell_pop_type is CellPop.gran:
            pop = list(map(float, list(np.array(config.cell_pop[cell_pop_type.value])[indexes])))

    return pop
