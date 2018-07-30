from config.config import *
from bop.approach.top.manova.top import save_top_manova


def top_proc(config, attributes_types, attribute_target):
    if config.approach_method is Method.manova:
        save_top_manova(config, attributes_types, attribute_target)


db = DataBaseType.GSE52588
dt = DataType.bop
approach = Approach.top
validation = Validation.simple
scenario = Scenario.approach
approach_method = Method.manova
validation_method = Method.linreg_mult
gender = Gender.any
disease = Disease.healthy
approach_gd = GeneDataType.mean
validation_gd = GeneDataType.mean
geo = GeoType.any
dna_region = DNARegion.any
cpg_class = ClassType.class_a
attributes_types = [Attribute.age]
attribute_target = Attribute.age



config = Config(
    db=db,
    dt=dt,
    approach=approach,
    validation=validation,
    scenario=scenario,
    approach_method=approach_method,
    validation_method=validation_method,
    gender=gender,
    disease=disease,
    approach_gd=approach_gd,
    validation_gd=validation_gd,
    geo=geo,
    dna_region=dna_region,
    cpg_class=cpg_class,
)

top_proc(config, attributes_types, attribute_target)
