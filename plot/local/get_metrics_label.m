function metrics_label = get_metrics_label(config)
metrics_label = '';
if strcmp(config.method, 'linreg')
    metrics_label = 'r';
elseif strcmp(config.method, 'manova')
    metrics_label = '$-\log_{10}p\_value$';
elseif strcmp(config.method, 'linreg_variance')
    metrics_label = 'r variance';
end
end