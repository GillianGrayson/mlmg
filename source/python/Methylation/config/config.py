from config.types.common import *

from config.types.annotations.annotation import *
from config.types.annotations.bop import *
from config.types.annotations.cpg import *
from config.types.annotations.gene import *
from config.types.annotations.common import *

from config.types.experiments.scenario import *
from config.types.experiments.approach import *
from config.types.experiments.method import *

from config.types.attributes.common import *
from config.types.attributes.attribute import *
from config.types.attributes.common import *
from config.types.attributes.cell_pop import *

from infrastructure.load.annotations import *
from infrastructure.load.attributes import *
from infrastructure.load.cell_pop import *
from infrastructure.load.excluded import *
from infrastructure.load.indexes import get_indexes
import socket
import getpass


class Config:

    def __init__(self,
                 read_only=False,

                 data_base=DataBase.GSE40279,
                 data_type=DataType.gene,

                 cross_reactive=CrossReactiveType.cross_reactive_included,
                 snp=SNPType.snp_included,
                 chromosome_type=ChromosomeType.non_gender,

                 class_type=ClassType.class_ab,

                 geo_type=GeoType.islands_shores,
                 gene_data_type=GeneDataType.mean,

                 dna_region=DNARegionType.any,

                 scenario=Scenario.approach,
                 approach=Approach.top,
                 method=Method.match,

                 disease=Disease.any,
                 gender=Gender.any,

                 attributes_types=None,
                 attribute_target=None,
                 method_params=None,

                 is_clustering=False
                 ):


        # Level 0: data_base
        self.data_base = data_base

        # Level 1: data_type
        self.data_type = data_type

        # Level 2: Annotations
        self.cross_reactive = cross_reactive
        self.snp = snp
        self.chromosome_type = chromosome_type
        # BOP
        self.class_type = class_type
        # GENE
        self.geo_type = geo_type
        self.gene_data_type = gene_data_type
        # CPG
        self.dna_region = dna_region

        # Level 3: Solution
        self.scenario = scenario
        self.approach = approach
        self.method = method

        # Level 4: Attributes
        self.disease = disease
        self.gender = gender

        self.attributes_types = attributes_types
        self.attribute_target = attribute_target

        self.method_params = method_params

        # Aux
        self.read_only = read_only

        # Clustering flag for methods
        self.is_clustering = is_clustering

        # Checking for GENE
        if self.data_type is DataType.gene:
            self.dna_region = DNARegionType.genic

        # Checking for GENE
        if self.data_type is DataType.cpg:
            self.geo_type = GeoType.any

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
            self.excluded = load_excluded(self)
            self.indexes = get_indexes(self)

        # Aux data
        self.print_rate = 10000
        self.num_skip_lines = 1
        self.miss_tag = 'NA'
        self.test_part = 0.25
        self.shift = 5
