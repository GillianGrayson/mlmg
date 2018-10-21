function metrics_id = get_metrics_id(config)
metrics_id = 1;

if strcmp(config.method, 'linreg')
    if config.metrics_rank == 1
        metrics_id = 1;
    elseif config.metrics_rank == 2
        metrics_id = 2;
     elseif config.metrics_rank == 3
        metrics_id = 3;
    end
    
elseif strcmp(config.method, 'manova')
    metrics_id = 1;
    
elseif strcmp(config.method, 'linreg_variance')
    if config.metrics_rank == 1
        metrics_id = 6;
    elseif config.metrics_rank == 2
        metrics_id = 1;
    end
    
elseif strcmp(config.method, 'moment')
    if config.metrics_rank == 1
        metrics_id = 1;
    elseif config.metrics_rank == 2
        metrics_id = 2;
    end    

elseif strcmp(config.method, 'linreg_ols')
    metrics_id = 1;

end
end