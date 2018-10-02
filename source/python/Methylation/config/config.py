from infrastructure.load.annotations import *
from infrastructure.load.attributes import *
from infrastructure.load.cell_pop import *
from infrastructure.load.indexes import get_indexes
import socket
import getpass


class Config:

    def __init__(self,
                 read_only=False,
                 data_base=DataBaseType.GSE40279,
                 data_type=DataType.gene,

                 chromosome_type=ChromosomeTypes.non_gender,

                 class_type=ClassType.class_ab,

                 gene_data_type=GeneDataType.mean,
                 geo_type=GeoType.islands_shores,

                 dna_region=DNARegion.any,

                 disease=Disease.any,
                 gender=Gender.any,
                 scenario=Scenario.approach,
                 approach=Approach.top,
                 method=Method.match,
                 ):
        # Global
        self.read_only = read_only
        self.data_base = data_base
        self.data_type = data_type

        # Common
        self.chromosome_type = chromosome_type

        # BOP
        self.class_type = class_type

        # GENE
        self.gene_data_type = gene_data_type
        self.geo_type = geo_type

        # CPG
        self.dna_region = dna_region

        # Experiment
        self.disease = disease
        self.gender = gender
        self.scenario = scenario
        self.approach = approach
        self.method = method

        # FS type
        host_name = socket.gethostname()
        self.fs = FSType.local_big
        if host_name == 'MSI':
            self.fs = FSType.local_msi
        elif host_name == 'DESKTOP-K9VO2TI':
            self.fs = FSType.local_big
        elif host_name == 'DESKTOP-4BEQ7MS':
            self.fs = FSType.local_ab
        elif host_name == 'master' or host_name[0:4] == 'node':
            user = getpass.getuser()
            if user == 'yusipov_i':
                self.fs = FSType.unn
            elif user == 'kalyakulina_a':
                self.fs = FSType.unn_ab
            else:
                self.fs = FSType.unn

        # Core data
        if not read_only:
            self.annotations = load_annotations(self)
            self.attributes = load_attributes(self)
            self.cell_pop = load_cell_pop(self)
            self.indexes = get_indexes(self)

        # Aux data
        self.print_rate = 10000
        self.num_skip_lines = 1
        self.miss_tag = 'NULL'
        self.test_part = 0.25
        self.shift = 5
