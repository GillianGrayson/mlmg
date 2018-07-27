from config.config import *
from annotation.bop import *
from infrastructure.load.cpg_data import load_cpg_data
from statsmodels.multivariate.manova import MANOVA

window = 3

db = DataBaseType.GSE40279
dt = DataType.gene
approach = Approach.top
validation = Validation.simple
scenario = Scenario.approach
approach_method = Method.manova
validation_method = Method.linreg_mult
gt = Gender.any
approach_gd = GeneDataType.mean
validation_gd = GeneDataType.mean
geo = GeoType.any
dna_region = DNARegion.any
cpg_class = ClassType.class_a

config = Config(
    db=db,
    dt=dt,
    approach=approach,
    validation=validation,
    scenario=scenario,
    approach_method=approach_method,
    validation_method=validation_method,
    gt=gt,
    approach_gd=approach_gd,
    validation_gd=validation_gd,
    geo=geo,
    dna_region=dna_region,
    cpg_class=cpg_class,
)

dict_bop_cpgs = get_dict_bop_cpgs(config)
cpgs, betas = load_cpg_data(config)
exog = np.array(get_main_attributes(config))

bops_passed = []
bops_pvals = []
for bop in dict_bop_cpgs:
    curr_cpgs = dict_bop_cpgs.get(bop)
    if len(curr_cpgs) > 2:
        for win_id in range(0, len(curr_cpgs) - 2):
            endog = []
            for cpg_id in range(0, window):
                cpg = curr_cpgs[win_id + cpg_id]
                beta = betas[cpgs.index(cpg)]
                endog.append(beta)
            endog = np.array(endog).T

            manova = MANOVA(endog, exog)
            mv_test_res =manova.mv_test()
            pfff = mv_test_res.results['x0']['stat'].values[0, 4]
            ololo = 1


a = 0
