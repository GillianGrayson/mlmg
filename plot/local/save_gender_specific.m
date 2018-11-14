function save_gender_specific(config, save_config)

suffix = sprintf('method(%s)', ...
    config.method);
path = sprintf('%s/data/%s', ...
    config.up, ...
    get_result_path(save_config));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

order = config.order;
names = config.names;
metrics_diff = config.metrics_diff;
metrics_diff_labels = config.metrics_diff_labels;

d = vertcat("names", string(names(order)));
for metrics_id = 1:size(metrics_diff, 2)
    d = horzcat(d, vertcat(metrics_diff_labels(metrics_id), string(metrics_diff(order, metrics_id))));
end

xlswrite(fn, d);

end