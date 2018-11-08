function gender_specific(config, save_config)

if strcmp(config.gender, 'any')
    [config.names, config.data] = get_data(config);
    save_data(config, save_config);
else
    [config.names, config.data_1, config.data_2] = get_gender_specific_data(config);
    [config.metrics_1, config.metrics_2] = get_gender_specific_metrics(config);
    [config.metrics_diff, config.metrics_diff_labels] = get_gender_specific_metrics_diff(config);
    
    num_names = size(config.names, 1)
    
    config.order = get_gender_specific_order(config);
    
    save_gender_specific(config, save_config);
    
    metrics_diff_srt = config.metrics_diff(config.order, config.metrics_diff_id);
    names_srt = config.names(config.order);
    metrics_1_srt = config.metrics_1(config.order);
    metrics_2_srt = config.metrics_2(config.order);
    
    plot_data.num_bins = 100;
    plot_data.num_rare = floor(config.part * size(names_srt, 1));
    plot_data.names = names_srt;
    plot_data.metrics_1 = metrics_1_srt;
    plot_data.metrics_2 = metrics_2_srt;
    plot_data.metrics_diff = metrics_diff_srt;
    plot_data.metrics_1_label = sprintf('%s F', get_metrics_label(config));
    plot_data.metrics_2_label = sprintf('%s M', get_metrics_label(config));
    plot_data.save_path = sprintf('%s/%s', ...
        save_config.up, ...
        get_result_path(save_config));
    plot_data.fig_name = 'gender';
    plot_data.suffix = sprintf('gender_method(%s)_rank(%d)_plot(%d)_part(%0.4f)', ...
        config.method, ...
        config.metrics_rank, ...
        config.plot_method, ...
        config.part);
    
    mkdir(plot_data.save_path);
    
    if strcmp(config.method, 'linreg_ols') || strcmp(config.method, 'linreg_ols_wo_outliers')
        plot_metrics_diff(plot_data);
    else
        plot_butterfly(plot_data);
        plot_metrics_diff(plot_data);
    end
end

end