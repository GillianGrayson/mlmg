function metrics_id = get_metrics_id(config, rank)
metrics_id = 1;

if strcmp(config.method, 'linreg')
    if rank == 1
        metrics_id = 1;
    elseif rank == 2
        metrics_id = 3;
    end
    
elseif strcmp(config.method, 'manova')
    metrics_id = 1;
    
elseif strcmp(config.method, 'linreg_variance')
    if rank == 1
        metrics_id = 5;
    elseif rank == 2
        metrics_id = 1;
    end
    
end
end