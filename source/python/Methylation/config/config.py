from config.types import *
from infrastructure.load.annotations import *
from infrastructure.load.attributes import *
from infrastructure.load.indexes import get_indexes
import socket


class Config:

    def __init__(self,
                 db=DataBaseType.GSE40279,
                 dt=DataType.gene,
                 approach=Approach.top,
                 validation=Validation.simple,
                 scenario=Scenario.approach,
                 approach_method=Method.linreg,
                 validation_method=Method.linreg_mult,
                 gender=Gender.any,
                 disease=Disease.any,
                 approach_gd=GeneDataType.mean,
                 validation_gd=GeneDataType.mean,
                 geo=GeoType.any,
                 dna_region=DNARegion.any,
                 cpg_class=ClassType.any,
                 ):
        # Config data
        self.db = db
        self.dt = dt
        self.approach = approach
        self.validation = validation
        self.scenario = scenario
        self.approach_method = approach_method
        self.validation_method = validation_method
        self.gender = gender
        self.disease = disease
        self.approach_gd = approach_gd
        self.validation_gd = validation_gd
        self.geo_type = geo
        self.class_type = cpg_class
        self.dna_region = dna_region

        # FS type
        host_name = socket.gethostname()
        self.fs = FSType.local_big
        if host_name == 'MSI':
            self.fs = FSType.local_msi
        elif host_name == 'DESKTOP-K9VO2TI':
            self.fs = FSType.local_big
        elif host_name == 'DESKTOP-4BEQ7MS':
            self.fs = FSType.local_ab

        # Core data
        self.annotations = load_annotations(self)
        self.attributes = load_attributes(self)
        self.indexes = get_indexes(self)

        # Aux data
        self.print_rate = 10000
        self.num_skip_lines = 1
        self.miss_tag = 'NULL'
        self.test_part = 0.25
        self.shift = 5
