function [passed_names, metrics_labels, metrics_data] = lvl_1_condition(config)
[names, data_1, data_2] = get_gender_specific_data(config);

if strcmp(config.data_type, 'linreg_ols')
    
    if experiment == 1
        
        sigma = 3;
        
        slopes_1 = data_1(:, 3);
        slopes_1_std = data_1(:, 5);
        slopes_2 = data_2(:, 3);
        slopes_2_std = data_2(:, 5);
        
        passed_names = [];
        for id = 1:size(names)
            if slopes_1(id) < sigma * slopes_1_std(id) && slopes_2(id) < sigma * slopes_2_std(id)
                passed_names = vertcat(passed_names, names(id));
            end
        end
        
    elseif config.experiment == 2
        
        slope_lim = 0.002;
        
        slopes_1 = data_1(:, 3);
        slopes_2 = data_2(:, 3);
        
        metrics_labels = ["slope_f"; "slope_m"];
        
        passed_names = [];
        metrics_data = []
        for id = 1:size(names)
            if abs(slopes_1(id)) > slope_lim || abs(slopes_2(id)) > slope_lim
                passed_names = vertcat(passed_names, names(id));
                metrics_data = vertcat(metrics_data, [slopes_1(id), slopes_2(id)])
            end
        end
        
    end
    
end

end