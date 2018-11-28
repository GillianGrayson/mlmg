function [passed_names, metrics_labels, metrics_map] = lvl_2_condition(method, config_lvl_2)

suffix = sprintf('method(%s)', ...
    method);
path = sprintf('%s/data/%s', ...
    config_lvl_2.up, ...
    get_result_path(config_lvl_2));
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

[num,txt,raw] = xlsread(fn);


if config_lvl_2.experiment == 1
    
    if strcmp(method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection = cell2mat(raw(2:end, 3));
        metrics_labels = [raw(1, 3)];
        passed_names = [];
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = area_intersection(id);
            end
        end
        
    end
    
elseif config_lvl_2.experiment == 2
    
    if strcmp(method, 'linreg_ols')
        
        names = raw(2:end, 1);
        slope_intersection = cell2mat(raw(2:end, 5));
        metrics_labels = [raw(1, 5)];
        passed_names = [];
        metrics_map = containers.Map(); 
        for id = 1:size(names)
            if slope_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = slope_intersection(id);
            end
        end
        
    end
    
elseif config_lvl_2.experiment == 3
    
    if strcmp(method, 'linreg_ols')
        
        names = raw(2:end, 1);
        slope_intersection = cell2mat(raw(2:end, 5));
        metrics_labels = [raw(1, 5)];
        passed_names = [];
        metrics_map = containers.Map(); 
        for id = 1:size(names)
            if slope_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = slope_intersection(id);
            end
        end
        
    end

elseif config_lvl_2.experiment == 4
    
    if strcmp(method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection = cell2mat(raw(2:end, 3));
        metrics_labels = [raw(1, 3)];
        passed_names = [];
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = area_intersection(id);
            end
        end
        
    end
    
end
    
end

