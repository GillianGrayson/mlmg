from method.linreg_mult.routines import *
from infrastructure.load.attributes import get_attributes
from infrastructure.load.top import load_top_gene_data
from infrastructure.file_system import *
from infrastructure.save.features import save_features

def save_simple_linreg_mult(config, num_bootstrap_runs=500, num_top=100):
    attributes = get_attributes(config)
    config.scenario = Scenario.approach
    gene_names, gene_vals = load_top_gene_data(config, num_top)
    config.scenario = Scenario.validation

    counts, R2s = R2_from_count(gene_vals, attributes)
    fn = 'R2s.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [counts, R2s])

    test_size = int(len(attributes) * config.test_part)
    train_size = len(attributes) - test_size
    metrics_names, metrics_vals = validation_metrics(gene_vals, attributes, test_size, train_size, num_bootstrap_runs)
    fn = 'metrics.txt'
    fn = get_result_path(config, fn)
    save_features(fn, [metrics_names, metrics_vals])

    print(linreg_mult_with_const(attributes, gene_vals).summary())

