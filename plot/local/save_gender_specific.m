function save_gender_specific(config_lvl_1, config_lvl_2)

suffix = sprintf('method(%s%s)', ...
    config_lvl_1.method, ...
    config_lvl_1.suffix);
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

order = config_lvl_1.order;
names = config_lvl_1.names;
metrics_diff = config_lvl_1.metrics_diff;
metrics_diff_labels = config_lvl_1.metrics_diff_labels;

d = vertcat("names", string(names(order)));
for metrics_id = 1:size(metrics_diff, 2)
    d = horzcat(d, vertcat(metrics_diff_labels(metrics_id), string(metrics_diff(order, metrics_id))));
end

xlswrite(fn, d);

end