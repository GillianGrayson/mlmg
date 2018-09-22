from config.config import *
from bop.approach.top.manova.top import save_top_manova


def top_proc(config, attributes_types, attribute_target):
    if config.approach_method is Method.manova:
        save_top_manova(config, attributes_types, attribute_target)


db = DataBaseType.GSE52588_TEST
dt = DataType.bop
scenario = Scenario.approach
approach = Approach.top
approach_method = Method.manova
gender = Gender.any
disease = Disease.any
cpg_classes = (ClassType.class_a)
attributes_types = [Attribute.age,
                    Attribute.gender,
                    CellPop.cd4_t,
                    CellPop.nk,
                    CellPop.b_cell,
                    CellPop.mono,
                    CellPop.gran,
                    CellPop.cd8_t,
                    Attribute.batch]
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
