function gender_specific(config_lvl_1, config_lvl_2)

if strcmp(config_lvl_1.gender, 'any')
    [config_lvl_1.names, config_lvl_1.data] = get_data(config_lvl_1);
    save_data(config_lvl_1, config_lvl_2);
else
    [config_lvl_1.names, config_lvl_1.data_1, config_lvl_1.data_2] = get_gender_specific_data(config_lvl_1);
    [config_lvl_1.metrics_1, config_lvl_1.metrics_2] = get_gender_specific_metrics(config_lvl_1);
    [config_lvl_1.metrics_diff, config_lvl_1.metrics_diff_labels] = get_gender_specific_metrics_diff(config_lvl_1);
    
    num_names = size(config_lvl_1.names, 1)
    
    config_lvl_1.order = get_gender_specific_order(config_lvl_1);
    
    save_gender_specific(config_lvl_1, config_lvl_2);
    
    metrics_diff_srt = config_lvl_1.metrics_diff(config_lvl_1.order, config_lvl_1.metrics_diff_id);
    names_srt = config_lvl_1.names(config_lvl_1.order);
    metrics_1_srt = config_lvl_1.metrics_1(config_lvl_1.order);
    metrics_2_srt = config_lvl_1.metrics_2(config_lvl_1.order);
    
    plot_data.num_bins = 100;
    plot_data.num_rare = floor(config_lvl_1.part * size(names_srt, 1));
    plot_data.names = names_srt;
    plot_data.metrics_1 = metrics_1_srt;
    plot_data.metrics_2 = metrics_2_srt;
    plot_data.metrics_diff = metrics_diff_srt;
    plot_data.metrics_1_label = sprintf('%s F', get_metrics_label(config_lvl_1));
    plot_data.metrics_2_label = sprintf('%s M', get_metrics_label(config_lvl_1));
    plot_data.save_path = sprintf('%s/%s', ...
        config_lvl_2.up, ...
        get_result_path(config_lvl_2));
    plot_data.fig_name = 'gender';
    plot_data.suffix = sprintf('gender_method(%s)_rank(%d)_plot(%d)_part(%0.4f)', ...
        config_lvl_1.method, ...
        config_lvl_1.metrics_rank, ...
        config_lvl_1.plot_method, ...
        config_lvl_1.part);
    
    mkdir(plot_data.save_path);
    
    if strcmp(config_lvl_1.method, 'linreg_ols')
        plot_metrics_diff(plot_data);
    else
        plot_butterfly(plot_data);
        plot_metrics_diff(plot_data);
    end
end

end