function metrics_label = get_metrics_label(config)
metrics_label = '';

if strcmp(config.method, 'linreg')
    if config.metrics_rank == 1
        metrics_label = 'r';
    elseif config.metrics_rank == 2
        metrics_label = 'slope';
    end
    
elseif strcmp(config.method, 'manova')
    metrics_label = '$-\log_{10}p\_value$';
    
elseif strcmp(config.method, 'linreg_variance')
    if config.metrics_rank == 1
        metrics_label = 'r variance';
    elseif config.metrics_rank == 2
        metrics_label = 'r';
    end

elseif strcmp(config.method, 'moment')
    if config.metrics_rank == 1
        metrics_label = 'mean';
    elseif config.metrics_rank == 2
        metrics_label = 'std';
    end
    
    
end
end
