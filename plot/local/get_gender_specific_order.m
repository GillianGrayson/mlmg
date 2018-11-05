function order = get_gender_specific_order(config, save_config)
    suffix = sprintf('method(%s)', ...
    config.method);
    path = sprintf('%s/data/%s', ...
        config.up, ...
        get_result_path(save_config));
    mkdir(path)
    fn = sprintf('%s/%s.xlsx', ...
        path, ...
        suffix);

if strcmp(config.method, 'linreg_ols')
    
    sigma = 3;
    
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
    areas_normed = zeros(size(names, 1), 1);
    for id = 1 : size(names, 1)
        
        intercept_minus_1 = intercepts_1(id) - sigma * intercepts_std_1(id);
        slope_minus_1 = slopes_1(id) - sigma * slopes_std_1(id);
        intercept_plus_1 = intercepts_1(id) + sigma * intercepts_std_1(id);
        slope_plus_1 = slopes_1(id) + sigma * slopes_std_1(id);
        
        intercept_up_1 = intercepts_1(id) + ((slope_plus_1 * x(2) + intercept_plus_1) - (slopes_1(id) * x(2) + intercepts_1(id)));
        intercept_down_1 = intercepts_1(id) + ((slope_minus_1 * x(2) + intercept_minus_1) - (slopes_1(id) * x(2) + intercepts_1(id)));
        
        y1 = [slopes_1(id) * x(1) + intercept_down_1, ...
            slopes_1(id) * x(2) + intercept_down_1, ...
            slopes_1(id) * x(3) + intercept_up_1, ...
            slopes_1(id) * x(4) + intercept_up_1];
        
        intercept_minus_2 = intercepts_2(id) - sigma * intercepts_std_2(id);
        slope_minus_2 = slopes_2(id) - sigma * slopes_std_2(id);
        intercept_plus_2 = intercepts_2(id) + sigma * intercepts_std_2(id);
        slope_plus_2 = slopes_2(id) + sigma * slopes_std_2(id);
        
        intercept_up_2 = intercepts_2(id) + ((slope_plus_2 * x(2) + intercept_plus_2) - (slopes_2(id) * x(2) + intercepts_2(id)));
        intercept_down_2 = intercepts_2(id) + ((slope_minus_2 * x(2) + intercept_minus_2) - (slopes_2(id) * x(2) + intercepts_2(id)));
        
        y2 = [slopes_2(id) * x(1) + intercept_down_2, ...
            slopes_2(id) * x(2) + intercept_down_2, ...
            slopes_2(id) * x(3) + intercept_up_2, ...
            slopes_2(id) * x(4) + intercept_up_2];
        
        pgon_1 = polyshape(x, y1);
        area_pgon_1 = polyarea(x, y1);
        pgon_2 = polyshape(x, y2);
        area_pgon_2 = polyarea(x, y2);
        
        pgon_intersect = intersect(pgon_1, pgon_2);
        
        areas(id) = polyarea(pgon_intersect.Vertices(:, 1), pgon_intersect.Vertices(:, 2));
        areas_normed(id) = areas(id) / (area_pgon_1 + area_pgon_2 - areas(id));
    end
    
    [tmp, order] = sort(areas, 'ascend');
    
    names = vertcat('names', names(order));
    areas = vertcat('areas', string(areas(order)));
    areas_normed = vertcat('areas_normed', string(areas_normed(order)));
    
    xlswrite(fn, [names, areas, areas_normed]);
    
else
    
    [tmp, order] = sort(abs(config.diff_metrics), 'descend');
    
end
end