function passed_names = lvl_2_condition(config_lvl_1, config_lvl_2)

suffix = sprintf('method(%s)', ...
    config_lvl_1.method);
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s_wo_cross_reactive.xlsx', ...
    path, ...
    suffix);

[num,txt,raw] = xlsread(fn);

names = raw(2:end, 1);
areas = cell2mat(raw(2:end, 3));

passed_names = [];
for id = 1:size(names)
    if areas(id) < 0.5
        passed_names = vertcat(passed_names, names(id));
    end
end

end

