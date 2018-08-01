from config.config import *
from bop.approach.top.manova.top import save_top_manova


def top_proc(config, attributes_types, attribute_target):
    if config.approach_method is Method.manova:
        save_top_manova(config, attributes_types, attribute_target)


db = DataBaseType.GSE40279
dt = DataType.bop
scenario = Scenario.approach
approach = Approach.top
approach_method = Method.manova
gender = Gender.M
disease = Disease.any
cpg_class = ClassType.class_a
attributes_types = [Attribute.age]
attribute_target = Attribute.age

config = Config(
    db=db,
    dt=dt,
    scenario=scenario,
    approach=approach,
    approach_method=approach_method,
    gender=gender,
    disease=disease,
    cpg_class=cpg_class,
)

top_proc(config, attributes_types, attribute_target)
