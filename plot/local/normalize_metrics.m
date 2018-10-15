function metrics_diff = normalize_metrics(target_metrics_diff, config)
metrics_diff = target_metrics_diff;
if strcmp(config.method, 'linreg')
    max_metrics_diff = max(abs(max(metrics_diff)), abs(min(metrics_diff)));
    min_metrics_diff = -max_metrics_diff;
    for gene_id = 1:size(metrics_diff, 1)
        metrics_diff(gene_id) = (metrics_diff(gene_id) - min_metrics_diff) / (max_metrics_diff - min_metrics_diff) * 2 - 1;
    end
end
end
