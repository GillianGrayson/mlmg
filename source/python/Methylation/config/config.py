from config.types import *
from infrastructure.load.annotations import *
from infrastructure.load.attributes import *
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
                 approach_gd=GeneDataType.mean,
                 validation_gd=GeneDataType.mean,
                 geo=GeoType.any,
                 dna_region=DNARegion.genic,
                 cpg_class=ClassType.any,
                 ):
        # Config data
        self.db = db
        self.dt = dt
        self.approach = approach
        self.validation = validation
        self.scenario=scenario
        self.approach_method = approach_method
        self.validation_method = validation_method
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

        # Core data
        self.annotations = load_annotations(self)
        self.attributes = load_attributes(self)

        # Aux data
        self.print_rate = 10000
        self.num_skip_lines = 0
        self.attribute_fn = ''
        self.miss_tag = ''
        self.train_size = 0
        self.test_size = 0
        if self.db is DataBaseType.GSE40279:
            self.num_skip_lines = 1
            self.attribute_fn = 'attributes.txt'
            self.miss_tag = 'NULL'
            self.train_size = 482
            self.test_size = 174
        elif self.db is DataBaseType.GSE52588:
            self.num_skip_lines = 87
            self.attribute_fn = 'attribute.txt'
            self.miss_tag = 'NULL'
            self.train_size = 0
            self.test_size = 0
