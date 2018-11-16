function [passed_names, metrics_labels, metrics_data] = lvl_2_condition(config_lvl_1, config_lvl_2)

suffix = sprintf('method(%s)', ...
    config_lvl_1.method);
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

[num,txt,raw] = xlsread(fn);

if strcmp(config.data_type, 'linreg_ols')
    
     if experiment == 1

        names = raw(2:end, 1);
        areas = cell2mat(raw(2:end, 3));
        metrics_labels = [raw(1, 3)];

        passed_names = [];
        metrics_data = [];
        for id = 1:size(names)
            if areas(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_data = vertcat(metrics_data, areas(id))
            end
        end
        
     elseif config.experiment == 2
         
        names = raw(2:end, 1);
        slope_intersection = cell2mat(raw(2:end, 5));
        metrics_labels = [raw(1, 5)];
        
        passed_names = [];
        metrics_data = [];
        for id = 1:size(names)
            if slope_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_data = vertcat(metrics_data, slope_intersection(id))
            end
        end

     end
        
end
end

