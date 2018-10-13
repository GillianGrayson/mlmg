function metrics_label = get_metrics_label(config)
metrics_label = '';
if strcmp(config.method, 'linreg')
    metrics_label = 'r';
elseif strcmp(config.method, 'manova')
    metrics_label = '-log_{10}p_value';
end
end