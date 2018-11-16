clear all;

config.experiment = 2;

% ======== config ========
config_lvl_1.data_base = 'GSE87571';
config_lvl_1.data_type = 'cpg_data';

config_lvl_1.cross_reactive = 'cross_reactive_excluded';
config_lvl_1.snp = 'snp_included';
config_lvl_1.chromosome_type = 'non_gender';

config_lvl_1.dna_region = 'genic';

config_lvl_1.info_type = 'result';

config_lvl_1.scenario = 'approach';
config_lvl_1.approach = 'top';
config_lvl_1.method = 'linreg_ols';

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

config_lvl_2.disease = config_lvl_1.disease;
config_lvl_2.gender = 'versus';

config_lvl_2.is_clustering = config_lvl_1.is_clustering;

config_lvl_2.up = config_lvl_1.up;

[lvl_1_names, lvl_1_metrics_labels, lvl_1_metrics_data] = lvl_1_condition(config_lvl_1);
[lvl_2_names, lvl_2_metrics_labels, lvl_2_metrics_data] = lvl_2_condition(config_lvl_1, config_lvl_2);

intersect_names = intersect(lvl_1_names, lvl_2_names);

metrics_data = [];
metrics_labels = [];

for name_id = 1:size(intersect_names,1)
    name = intersect_names(name_id);
    
    
end






order = config.order;
names = config.names;
metrics_diff = config.metrics_diff;
metrics_diff_labels = config.metrics_diff_labels;

d = vertcat("names", string(names(order)));
for metrics_id = 1:size(metrics_diff, 2)
    d = horzcat(d, vertcat(metrics_diff_labels(metrics_id), string(metrics_diff(order, metrics_id))));
end

xlswrite(fn, d);






suffix = sprintf('lvl_3_experiment(%s)', config.experiment);
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

d = vertcat('names', intersect_names);

xlswrite(fn, d);