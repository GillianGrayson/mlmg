function order = get_gender_specific_order(config)
if strcmp(config.method, 'linreg_ols')
    [tmp, order] = sort(abs(config.diff_metrics), 'descend');
else
    [tmp, order] = sort(abs(config.diff_metrics), 'descend');
end
end