from config.config import *
from bop.approach.top.manova.top import save_top_manova


def top_proc(config, attributes_types, attribute_target):
    if config.approach_method is Method.manova:
        save_top_manova(config, attributes_types, attribute_target)


db = DataBaseType.GSE30870
dt = DataType.bop
scenario = Scenario.approach
approach = Approach.top
approach_method = Method.manova
gender = Gender.any
disease = Disease.any
cpg_classes = (ClassType.class_a, ClassType.class_b)
attributes_types = [Attribute.age, CellPop.plasma_blast, CellPop.cd8_p, CellPop.cd8_naive, CellPop.cd4_naive,
                    CellPop.cd8_t, CellPop.cd4_t, CellPop.nk, CellPop.b_cell, CellPop.mono, CellPop.gran]
attribute_target = Attribute.age

config = Config(
    db=db,
    dt=dt,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=gender,
    disease=disease,
    cpg_classes=cpg_classes,
)

top_proc(config, attributes_types, attribute_target)
