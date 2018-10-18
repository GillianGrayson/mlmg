function [metrics_1, metrics_2] = process_metrics_plane(target_metrics_1, target_metrics_2, config)
metrics_1 = target_metrics_1;
metrics_2 = target_metrics_2;
if strcmp(config.method, 'manova')
    if config.plot_method == 2
        p_value_lim = 1e-8;
        p_value_lim_log = -log10(p_value_lim);
        metrics_1 = [];
        metrics_2 = [];
        for id = 1:size(target_metrics_1, 1)
            if(target_metrics_1(id) < p_value_lim_log) || (target_metrics_2(id) < p_value_lim_log)
                metrics_1 = vertcat(metrics_1, target_metrics_1(id));
                metrics_2 = vertcat(metrics_2, target_metrics_2(id));
            end
        end
    end
elseif strcmp(config.method, 'linreg')
    if config.plot_method == 2
        p_value_lim = 1e-8;
        p_value_lim_log = -log10(p_value_lim);
        metrics_1 = [];
        metrics_2 = [];
        for id = 1:size(target_metrics_1, 1)
            if(target_metrics_1(id) < p_value_lim_log) || (target_metrics_2(id) < p_value_lim_log)
                metrics_1 = vertcat(metrics_1, target_metrics_1(id));
                metrics_2 = vertcat(metrics_2, target_metrics_2(id));
            end
        end
    end
end
end