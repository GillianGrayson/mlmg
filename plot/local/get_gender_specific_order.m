function order = get_gender_specific_order(config)
if strcmp(config.method, 'linreg_ols')
    
    config.names;
    config.f_metrics;
    config.m_metrics;
    config.diff_metrics;
    config.data_1;
    config.data_2;
    
    ages = get_ages(config);
    
    names = config.names;
    
    intercepts_1 = config.data_1(:, 2);
    slopes_1 = config.data_1(:, 3);
    intercepts_std_1 = config.data_1(:, 4);
    slopes_std_1 = config.data_1(:, 5);
    
    intercepts_2 = config.data_2(:, 2);
    slopes_2 = config.data_2(:, 3);
    intercepts_std_2 = config.data_2(:, 4);
    slopes_std_2 = config.data_2(:, 5);
    
    x = [min(ages), max(ages), max(ages), min(ages)];
    
    areas = zeros(size(names, 1), 1);
    for id = 1 : size(names, 1)
        intercept_down_1 = intercepts_1(id) - 3 * intercepts_std_1(id);
        slope_down_1 = slopes_1(id) - 3 * slopes_std_1(id);
        intercept_up_1 = intercepts_1(id) + 3 * intercepts_std_1(id);
        slope_up_1 = slopes_1(id) + 3 * slopes_std_1(id);
        y1 = [slope_down_1 * x(1) + intercept_down_1, ...
            slope_down_1 * x(2) + intercept_down_1, ...
            slope_up_1 * x(3) + intercept_up_1, ...
            slope_up_1 * x(4) + intercept_up_1];
        
        intercept_down_2 = intercepts_2(id) - 3 * intercepts_std_2(id);
        slope_down_2 = slopes_2(id) - 3 * slopes_std_2(id);
        intercept_up_2 = intercepts_2(id) + 3 * intercepts_std_2(id);
        slope_up_2 = slopes_2(id) + 3 * slopes_std_2(id);
        y2 = [slope_down_2 * x(1) + intercept_down_2, ...
            slope_down_2 * x(2) + intercept_down_2, ...
            slope_up_2 * x(3) + intercept_up_2, ...
            slope_up_2 * x(4) + intercept_up_2];
        
        pgon_1 = polyshape(x, y1);
        pgon_2 = polyshape(x, y2);
        
        pgon_intersect = intersect(pgon_1, pgon_2);
        
        
    end
    
    
    
    [tmp, order] = sort(abs(config.diff_metrics), 'descend');
    
else
    
    [tmp, order] = sort(abs(config.diff_metrics), 'descend');
    
end
end