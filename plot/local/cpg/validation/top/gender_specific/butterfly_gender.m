clear all;

% ======== params ========
config.metrics_rank = 3;
config.plot_method = 2;
config.part = 0.0005;

plot_data.num_bins = 200;

% ======== config ========
config.data_base = 'GSE40279';
config.data_type = 'cpg_data';

config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg';

config.disease = 'any';
config.gender = '';

config.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI')
    config.up = 'D:/YandexDisk/Work/mlmg';
else
    config.up = 'E:/YandexDisk/Work/mlmg';
end

% ======== save_config ========
save_config.data_base = config.data_base;
save_config.data_type = config.data_type;

save_config.chromosome_type = config.chromosome_type;

save_config.dna_region = config.dna_region;

save_config.info_type = 'result';

save_config.scenario = 'validation';
save_config.approach = 'top';
save_config.method = 'gender_specific';

save_config.disease = config.disease;
save_config.gender = 'versus';

save_config.up = 'C:/Users/user/Google Drive/mlmg/figures';
save_config.is_clustering = config.is_clustering;

% ======== processing ========
[names, data_1, data_2] = get_gender_specific_data(config);
config.names = names;
config.data_1 = data_1;
config.data_2 = data_2;
[metrics_1, metrics_2] = get_gender_specific_metrics(config);
config.metrics_1 = metrics_1;
config.metrics_2 = metrics_2;
diff_metrics = get_gender_specific_diff_metrics(config);
config.diff_metrics = diff_metrics;

order = get_gender_specific_order(config);

diff_metrics_srt = diff_metrics(order);
cpgs_srt = names(order);
metrics_1_srt = f_metrics(order);
metrics_2_srt = m_metrics(order);

plot_data.num_rare = floor(config.part * size(cpgs_srt, 1));
plot_data.names = cpgs_srt;
plot_data.metrics_1 = metrics_1_srt;
plot_data.metrics_2 = metrics_2_srt;
plot_data.metrics_diff = diff_metrics_srt;
plot_data.metrics_1_label = sprintf('%s F', get_metrics_label(config));
plot_data.metrics_2_label = sprintf('%s M', get_metrics_label(config));
plot_data.save_path = sprintf('%s/%s', ...
    save_config.up, ...
    get_result_path(save_config));
plot_data.fig_name = 'gender';
plot_data.suffix = sprintf('gender_method(%s)_rank(%d)_plot(%d)_part(%0.4f)', ...
    config.method, ...
    config.metrics_rank, ...
    config.plot_method, ...
    config.part);

mkdir(plot_data.save_path);

plot_butterfly(plot_data);
