function metrics_label = get_metrics_label(config, rank)
metrics_label = '';

if strcmp(config.method, 'linreg')
    if rank == 1
        metrics_label = 'r';
    elseif rank == 2
        metrics_label = 'slope';
    end
    
elseif strcmp(config.method, 'manova')
    metrics_label = '$-\log_{10}p\_value$';
    
elseif strcmp(config.method, 'linreg_variance')
    if rank == 1
        metrics_label = 'r variance';
    elseif rank == 2
        metrics_label = 'r';
    end
    
end
end
