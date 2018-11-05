function order = get_gender_specific_order(config)
    [tmp, order] = sort(config.diff_metrics(config.diff_metrics_id), config.diff_metrics_direction);
end