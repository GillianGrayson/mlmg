from enum import Enum


class FSType(Enum):
    local_big = 'E:/Work/mlmg/data'
    local_msi = 'D:/Work/mlmg/data'
    local_ab = 'D:/Aaron/Bio/mlmg/data'
    unn = '/common/home/yusipov_i/Work/mlmg/data'
    mpipks = '/data/biophys/yusipov/mlmg/data'
    unn_ab = '/common/home/kalyakulina_a/Work/mlmg/data'


class DataBaseType(Enum):
    GSE40279 = 'GSE40279'
    GSE52588 = 'GSE52588'
    GSE30870 = 'GSE30870'
    GSE61256 = 'GSE61256'
    GSE52588_TEST = 'GSE52588_TEST'


class DataType(Enum):
    cpg = 'cpg'
    gene = 'gene'
    bop = 'bop'


class Approach(Enum):
    top = 'top'
    clustering = 'clustering'
    bend = 'bend'
    statistics = 'statistics'


class Validation(Enum):
    simple = 'simple'


class Scenario(Enum):
    approach = 'approach'
    validation = 'validation'


class Gender(Enum):
    M = 'M'
    F = 'F'
    any = 'any'

class Disease(Enum):
    any = 'any'
    healthy = 'healthy'
    down_syndrome = 'down_syndrome'


class GeneDataType(Enum):
    mean = 'mean'
    std = 'std'
    mean_der = 'mean_der'
    mean_der_normed = 'mean_der_normed'
    from_cpg = 'from_cpg'
    from_bop = 'from_bop'


class GeoType(Enum):
    shores = 'shores'
    shores_s = 'shores_s'
    shores_n = 'shores_n'
    islands = 'islands'
    islands_shores = 'islands_shores'
    any = 'any'


class ClassType(Enum):
    class_a = 'ClassA'
    class_b = 'ClassB'
    class_c = 'ClassC'
    class_d = 'ClassD'


class DNARegion(Enum):
    genic = 'genic'
    non_genic = 'non_genic'
    any = 'any'


class Method(Enum):
    enet = 'enet'
    linreg = 'linreg'
    linreg_with_rejection = 'linreg_with_rejection'
    anova = 'anova'
    spearman = 'spearman'
    manova = 'manova'
    random_forest = 'random_forest'
    k_means = 'k_means'
    mean_shift = 'mean_shift'
    linreg_mult = 'linreg_mult'
    match = 'match'


class MANOVATest(Enum):
    wilks = 'wilks'
    pillai_bartlett = 'pillai_bartlett'
    lawley_hotelling = 'lawley_hotelling'
    roy = 'roy'


class Attribute(Enum):
    geo_accession = 'geo_accession'
    age = 'age'
    gender = 'gender'
    disease = 'disease'
    group = 'Group'
    batch = 'Batch'


class Annotation(Enum):
    cpg = 'ID_REF'
    chr = 'CHR'
    map_info = 'MAPINFO'
    gene = 'UCSC_REFGENE_NAME'
    class_type = 'Class'
    geo = 'RELATION_TO_UCSC_CPG_ISLAND'
    bop = 'BOP'


class CellPop(Enum):
    sample_id = 'SampleID'
    plasma_blast = 'PlasmaBlast'
    cd8_p = 'CD8pCD28nCD45RAn'
    cd8_naive = 'CD8naive'
    cd4_naive = 'CD4naive'
    cd8_t = 'CD8T'
    cd4_t = 'CD4T'
    nk = 'NK'
    b_cell = 'Bcell'
    mono = 'Mono'
    gran = 'Gran'

class CpGCondition(Enum):
    regular = 'regular'
    x = 'x'