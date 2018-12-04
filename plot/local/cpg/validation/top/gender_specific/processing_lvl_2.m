clear all;

% ======== params ========
config_lvl_1.metrics_rank = 1;
config_lvl_1.plot_method = 1;
config_lvl_1.metrics_diff_id = 2;
config_lvl_1.metrics_diff_direction = 'ascend';
config_lvl_1.part = 0.0005;

% ======== config ========
config_lvl_1.data_base = 'GSE87571';
config_lvl_1.data_type = 'cpg_data';

config_lvl_1.cross_reactive = 'cross_reactive_excluded';
config_lvl_1.snp = 'snp_excluded';
config_lvl_1.chromosome_type = 'non_gender';

config_lvl_1.dna_region = 'genic';

config_lvl_1.info_type = 'result';

config_lvl_1.scenario = 'approach';
config_lvl_1.approach = 'top';
config_lvl_1.method = 'linreg_ols';
config_lvl_1.suffix = '_outliers_limit(0.8)_outliers_sigma_(3.0)';

config_lvl_1.disease = 'any';
config_lvl_1.gender = 'versus';

config_lvl_1.is_clustering = 0;

config_lvl_1.up = get_up_data_path(); 

% ======== save_config ========
config_lvl_2.data_base = config_lvl_1.data_base;
config_lvl_2.data_type = config_lvl_1.data_type;

config_lvl_2.cross_reactive = config_lvl_1.cross_reactive;
config_lvl_2.snp = config_lvl_1.snp;
config_lvl_2.chromosome_type = config_lvl_1.chromosome_type;

config_lvl_2.dna_region = config_lvl_1.dna_region;

config_lvl_2.info_type = 'result';

config_lvl_2.scenario = 'validation';
config_lvl_2.approach = 'top';
config_lvl_2.method = 'gender_specific';
config_lvl_2.suffix = config_lvl_1.suffix;

config_lvl_2.disease = config_lvl_1.disease;
config_lvl_2.gender = 'versus';

config_lvl_2.is_clustering = config_lvl_1.is_clustering;

config_lvl_2.up = get_up_figures_path(); 

% ======== processing ========
gender_specific(config_lvl_1, config_lvl_2);