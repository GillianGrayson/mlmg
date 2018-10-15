function metrics_id = get_metrics_id(config)
metrics_id = 1;
if strcmp(config.method, 'linreg')
    metrics_id = 1;
elseif strcmp(config.method, 'manova')
    metrics_id = 1;
elseif strcmp(config.method, 'linreg_variance')
    metrics_id = 5;
end
end