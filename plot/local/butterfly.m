function butterfly(config, save_config)

[names, data_1, data_2] = get_gender_specific_data(config);
config.names = names;
config.data_1 = data_1;
config.data_2 = data_2;
[metrics_1, metrics_2] = get_gender_specific_metrics(config);
config.metrics_1 = metrics_1;
config.metrics_2 = metrics_2;
diff_metrics = get_gender_specific_diff_metrics(config);
config.diff_metrics = diff_metrics;

num_names = size(config.names, 1)

order = get_gender_specific_order(config);

diff_metrics_srt = diff_metrics(order);
names_srt = names(order);
metrics_1_srt = metrics_1(order);
metrics_2_srt = metrics_2(order);

plot_data.num_bins = 100;
plot_data.num_rare = floor(config.part * size(names_srt, 1));
plot_data.names = names_srt;
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

end