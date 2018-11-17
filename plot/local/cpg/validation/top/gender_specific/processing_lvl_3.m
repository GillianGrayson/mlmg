clear all;

config_lvl_1.experiment = 2;
config_lvl_2.experiment = 1;

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

[lvl_1_names, lvl_1_metrics_labels, lvl_1_metrics_map] = lvl_1_condition(config_lvl_1);
[lvl_2_names, lvl_2_metrics_labels, lvl_2_metrics_map] = lvl_2_condition(config_lvl_1, config_lvl_2);

intersect_names = intersect(lvl_1_names, lvl_2_names);

metrics_data = [];
metrics_labels = horzcat(lvl_1_metrics_labels, lvl_2_metrics_labels);

for name_id = 1:size(intersect_names,1)
    name = string(intersect_names(name_id));
    data = horzcat(lvl_1_metrics_map(name), lvl_2_metrics_map(name));
    metrics_data = vertcat(metrics_data, data);
end

suffix = sprintf('lvl_3_by_lvl_1(%d)_lvl_2(%d)', config_lvl_1.experiment, config_lvl_2.experiment);
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

d = vertcat("names", intersect_names);
for metrics_id = 1:size(metrics_data, 2)
    d = horzcat(d, vertcat(metrics_labels(metrics_id), string(metrics_data(:, metrics_id))));
end

xlswrite(fn, d);