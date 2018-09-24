from method.linreg_mult.routines import *
from infrastructure.load.top import load_top_gene_data
from infrastructure.save.features import save_features
from config.config import *

def save_error_from_age(config, num_top=100):
    attributes = get_attributes(config)
    config.scenario = Scenario.approach
    names, vals = load_top_gene_data(config, num_top)
    config.scenario = Scenario.validation

    X = vals
    y = attributes

    model = linreg_mult(y, X)

    ages = []
    maes = []
    str_list = []
    x_all = []
    y_all = []
    for age in range(0, 150):

        indexes = [i for i, x in enumerate(attributes) if x == age]

        if len(indexes) > 0:

            ages.append(age)

            X_test = np.array(vals).T[indexes].tolist()
            y_test_pred = model.get_prediction(X_test).predicted_mean

            curr_str = str(age)
            mae = 0
            for pred_age in y_test_pred:
                mae += abs(pred_age - age)
                curr_str += (' ' + str(format(pred_age, '0.8e')))
                x_all.append(age)
                y_all.append(pred_age - age)

            mae /= len(indexes)

            str_list.append(curr_str)
            maes.append(mae)

    fn = 'error_from_age.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [ages, maes])

    fn = 'errors.txt'
    fn = get_result_path(config, fn)
    np.savetxt(fn, str_list, fmt="%s")

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_all, y_all)
    print('slope: ' + str(slope))
    print('intercept: ' + str(intercept))
    print('r_value: ' + str(r_value))
    print('p_value: ' + str(p_value))
    print('std_err: ' + str(std_err))


num_top = 100
genders = [Gender.F, Gender.any, Gender.M]

for gender in genders:

    print('gender: ' + gender.value)

    config = Config(
        db=DataBaseType.GSE40279,
        dt=DataType.gene,
        scenario=Scenario.validation,
        approach=Approach.top,
        validation=Validation.simple,
        approach_method=Method.enet,
        validation_method=Method.linreg_mult,
        approach_gd=GeneDataType.mean,
        validation_gd=GeneDataType.mean,
        geo=GeoType.islands_shores,
        gender=gender
    )

    save_error_from_age(config, num_top=num_top)

