function [passed_names, metrics_labels, metrics_map] = lvl_2_condition(config_lvl_1, config_lvl_2)

suffix = sprintf('method(%s%s)', ...
    config_lvl_1.method, ...
    config_lvl_1.suffix);
path = sprintf('%s/data/%s', ...
    config_lvl_2.up, ...
    get_result_path(config_lvl_2));
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

[num,txt,raw] = xlsread(fn);

if config_lvl_2.experiment == 1
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection_rel = cell2mat(raw(2:end, 3));
        metrics_labels = [raw(1, 3)];
        passed_names = [];
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection_rel(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = area_intersection_rel(id);
            end
        end
        
    end
    
elseif config_lvl_2.experiment == 2
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        
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
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection_rel = cell2mat(raw(2:end, 3));
        slope_intersection = cell2mat(raw(2:end, 5));
        metrics_labels = [raw(1, 3), raw(1, 5)];
        passed_names = [];
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection_rel(id) < 0.5 && slope_intersection(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = [area_intersection_rel(id), slope_intersection(id)];
            end
        end
        
    end
    
elseif config_lvl_2.experiment == 4
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection_rel = cell2mat(raw(2:end, 3));
        metrics_labels = [raw(1, 3)];
        passed_names = [];
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection_rel(id) < 0.5
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = area_intersection_rel(id);
            end
        end
        
    end
    
elseif config_lvl_2.experiment == 5
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        
        names = raw(2:end, 1);
        area_intersection_rel = cell2mat(raw(2:end, 3));
        variance = cell2mat(raw(2:end, 4));
        metrics_labels = [raw(1, 3), raw(1, 4)];
        passed_names = strings(size(names, 1), 1);
        num_names = 1;
        metrics_map = containers.Map();
        for id = 1:size(names)
            if area_intersection_rel(id) < 0.5 && variance(id) > 3.0
                passed_names(num_names) = names(id);
                metrics_map(string(names(id))) = [area_intersection_rel(id), variance(id)];
                num_names = num_names + 1;
            end
        end
        num_names = num_names - 1;
        
        passed_names = passed_names(1:num_names, :);
        
    end
    
elseif config_lvl_2.experiment == 6
    
    if strcmp(config_lvl_1.method, 'linreg_variance_ols')
        names = raw(2:end, 1);
        slope_intersection = cell2mat(raw(2:end, 5));
        slope_intersection_var = cell2mat(raw(2:end, 9));
        metrics_labels = [raw(1, 5), raw(1, 9)];
        passed_names = strings(size(names, 1), 1);
        num_names = 1;
        metrics_map = containers.Map();
        for id = 1:size(names)
            passed_names(num_names) = names(id);
            metrics_map(string(names(id))) = [slope_intersection(id), slope_intersection_var(id)];
            num_names = num_names + 1;
        end
        num_names = num_names - 1;
        passed_names = passed_names(1:num_names, :);
    end
    
elseif config_lvl_2.experiment == 7
    
    if strcmp(config_lvl_1.method, 'linreg_variance_ols')
        names = raw(2:end, 1);
        variance = cell2mat(raw(2:end, 4));
        slope_intersection = cell2mat(raw(2:end, 5));
        slope_intersection_var = cell2mat(raw(2:end, 9));
        metrics_labels = [raw(1, 4), raw(1, 5), raw(1, 9)];
        passed_names = strings(size(names, 1), 1);
        num_names = 1;
        metrics_map = containers.Map();
        for id = 1:size(names)
            passed_names(num_names) = names(id);
            metrics_map(string(names(id))) = [variance(id), slope_intersection(id), slope_intersection_var(id)];
            num_names = num_names + 1;
        end
        num_names = num_names - 1;
        passed_names = passed_names(1:num_names, :);
    end

    
end

end

