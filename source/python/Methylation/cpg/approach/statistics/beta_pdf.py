from infrastructure.load.cpg_data import load_dict_cpg_data
from infrastructure.save.features import save_features
from config.config import *
import math


db = DataBase.GSE52588
dt = DataType.cpg
approach = Approach.statistics
scenario = Scenario.approach
genders = [Gender.any, Gender.M, Gender.F]
diseases = [Disease.any, Disease.down_syndrome, Disease.healthy]
geos = [GeoType.islands_shores, GeoType.islands, GeoType.any]
cpg_condition = CpGCondition.x

for disease in diseases:
    print('disease: ' + disease.value)
    for gender in genders:
        print('\tgender: ' + gender.value)
        for geo in geos:
            print('\t\tgeo: ' + geo.value)

            config = Config(db=db,
                            dt=dt,
                            approach=approach,
                            scenario=scenario,
                            gender=gender,
                            disease=disease,
                            geo=geo,
                            cpg_condition=cpg_condition)

            attributes = get_attributes(config)
            dict_cpg_data = load_dict_cpg_data(config)
            cpgs = list(dict_cpg_data.keys())
            vals = list(dict_cpg_data.values())

            num_int = 200
            int_begin = 0
            int_end = 1
            int_shift = (int_end - int_begin) / num_int
            ints = []
            pdf = np.zeros(num_int)
            for int_id  in range(0, num_int):
                ints.append(int_begin + int_id * int_shift + 0.5 * int_shift)

            for curr_cpg_vals in vals:
                for beta in curr_cpg_vals:
                    int_id = math.floor((beta - int_begin) * num_int / (int_end - int_begin + 1.0e-8))
                    pdf[int_id] += 1

            pdf = np.asarray(pdf)
            sum_pdf = np.sum(pdf)
            pdf = pdf / (sum_pdf * int_shift)
            print('pdf norm: ' + str(np.sum(pdf) * int_shift))

            fn = 'top.txt'
            fn = get_result_path(config, fn)
            save_features(fn, [ints, pdf])
            config.dt = DataType.cpg
