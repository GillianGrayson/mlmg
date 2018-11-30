clear all;

% ======== params ========
config.metrics_rank = 1;
config.plot_method = 1;
config.metrics_diff_id = 2;
config.metrics_diff_direction = 'ascend';
config.part = 0.0005;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'cpg_data';

config.cross_reactive = 'cross_reactive_excluded';
config.snp = 'snp_excluded';
config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_variance_ols';

config.disease = 'any';
config.gender = 'versus';

config.suffix = '';

config.is_clustering = 0;

config.up = get_up_data_path(); 

% ======== save_config ========
save_config.data_base = config.data_base;
save_config.data_type = config.data_type;

save_config.cross_reactive = config.cross_reactive;
save_config.snp = config.snp;
save_config.chromosome_type = config.chromosome_type;

save_config.dna_region = config.dna_region;

save_config.info_type = 'result';

save_config.scenario = 'validation';
save_config.approach = 'top';
save_config.method = 'gender_specific';

save_config.disease = config.disease;
save_config.gender = 'versus';

save_config.is_clustering = config.is_clustering;

save_config.up = get_up_figures_path(); 

% ======== processing ========
gender_specific(config, save_config);