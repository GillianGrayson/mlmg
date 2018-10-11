function metrics_id = get_metrics_id(config)
metrics_id = 1;
if strcmp(config.method, 'linreg')
    if config.is_clustering == 1
        metrics_id = 1;
    else
        metrics_id = 1;
    end
elseif strcmp(config.method, 'manova')
    if config.is_clustering == 1
        metrics_id = 1;
    else
        metrics_id = 1;
    end
end
end