function [passed_names, metrics_labels, metrics_map] = lvl_1_condition(config)
[names, data_1, data_2] = get_gender_specific_data(config);

if strcmp(config.method, 'linreg_ols')
    
    if config.experiment == 1
        
        sigma = 3;
        
        slopes_1 = data_1(:, 3);
        slopes_1_std = data_1(:, 5);
        slopes_2 = data_2(:, 3);
        slopes_2_std = data_2(:, 5);
        
        metrics_labels = ["slope_f", "slope_m"];
        
        passed_names = [];
        metrics_map = containers.Map(); 
        for id = 1:size(names)
            if slopes_1(id) < sigma * slopes_1_std(id) && slopes_2(id) < sigma * slopes_2_std(id)
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = [slopes_1(id), slopes_2(id)];
            end
        end
        
    elseif config.experiment == 2
        
        slope_lim = 0.002;
        
        slopes_1 = data_1(:, 3);
        slopes_2 = data_2(:, 3);
        slope_pvals_1 = data_1(:, 7);
        slope_pvals_2 = data_2(:, 7);
        
        metrics_labels = ["slope_f", "slope_m", "slope_pvals_f", "slope_pvals_m"];
        
        passed_names = [];
        metrics_map = containers.Map(); 
        for id = 1:size(names)
            if abs(slopes_1(id)) > slope_lim || abs(slopes_2(id)) > slope_lim
                passed_names = vertcat(passed_names, names(id));
                metrics_map(string(names(id))) = [slopes_1(id), slopes_2(id), slope_pvals_1(id), slope_pvals_2(id)];
            end
        end
        
    end
    
end

end