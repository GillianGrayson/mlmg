function metrics = process_metrics(target_metrics, config)
metrics = target_metrics;
if strcmp(config.method, 'manova')
    metrics = -log10(target_metrics);
end
end
