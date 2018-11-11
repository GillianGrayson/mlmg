function [metrics_diff, metrics_diff_labels] = get_gender_specific_metrics_diff(config)

if strcmp(config.method, 'linreg_ols') || strcmp(config.method, 'linreg_ols_wo_outliers')
    
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
    variance_diff = zeros(size(names, 1), 1);
    slope_intersection = zeros(size(names, 1), 1);
    
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
        
        variance_1 = y1(4) - y1(1);
        variance_2 = y2(4) - y2(1);
        variance_diff(id) = max(variance_1, variance_2) / min(variance_1, variance_2);
        
        pgon_slope_1_x = [slope_minus_1, slope_plus_1, slope_plus_1, slope_minus_1];
        pgon_slope_1_y = [0.0, 0.0, 1.0, 1.0];
        pgon_slope_1 = polyshape(pgon_slope_1_x, pgon_slope_1_y);
        
        pgon_slope_2_x = [slope_minus_2, slope_plus_2, slope_plus_2, slope_minus_2];
        pgon_slope_2_y = [0.0, 0.0, 1.0, 1.0];
        pgon_slope_2 = polyshape(pgon_slope_2_x, pgon_slope_2_y);
        
        pgon_slope_intersect = intersect(pgon_slope_1, pgon_slope_2);
        pgon_slope_union = union(pgon_slope_1, pgon_slope_2);
        
        area_intersection = polyarea(pgon_slope_intersect.Vertices(:, 1), pgon_slope_intersect.Vertices(:, 2));
        area_union = polyarea(pgon_slope_union.Vertices(:, 1), pgon_slope_union.Vertices(:, 2));
        
        if area_intersection < 1e-8
            slope_intersection(id) = 0.0;
        else
            slope_intersection(id) = area_intersection / area_union;
        end
        
    end
    
    metrics_diff = horzcat(areas, areas_normed, variance_diff, slope_intersection);
    metrics_diff_labels = ["areas", "areas_normed", "variance", "slope_intersection"];
    
else
    
    metrics_diff = zeros(size(config.names, 1), 1);
    for gene_id = 1:size(config.names, 1)
        metrics_diff(gene_id) = abs(config.metrics_1(gene_id) - config.metrics_2(gene_id));
    end
    metrics_diff_labels = ["metric"];

end