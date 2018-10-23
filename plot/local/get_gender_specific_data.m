function [names, data_1, data_2] = get_gender_specific_data(config)

config.gender = 'F';
fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));
raw_data_1 = importdata(fn, ' ');

config.gender = 'M';
fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));
raw_data_2 = importdata(fn, ' ');

names_1 = raw_data_1.textdata;
d_1_tmp = raw_data_1.data;

names_2 = raw_data_2.textdata;
d_2_tmp = raw_data_2.data;

all_names = names_1;
num_names = size(all_names, 1);

all_data_1 = d_1_tmp;
all_data_2 = zeros(num_names, 1);
for id_1 = 1:num_names
    name = string(all_names(id_1));
    id_2 = find(names_2==name);
    for m_id = 1:size(all_data_1, 2)
        all_data_2(id_1, m_id) = d_2_tmp(id_2, m_id);
    end
    
    if mod(id_1, 1000) == 0
        id_1 = id_1
    end
end

names = [];
data_1 = [];
data_2 = [];

if strcmp(config.method, 'linreg') && config.plot_method == 2
    
    p_value_lim = 1e-8;
    p_value_lim_log = -log10(p_value_lim);
    
    config_1 = config;
    config_1.metrics_rank = 2;
    
    metrics_id = 2;
    metrics_1 = all_data_1(:, metrics_id);
    metrics_2 = all_data_2(:, metrics_id);
    metrics_1 = process_metrics(metrics_1, config_1);
    metrics_2 = process_metrics(metrics_2, config_1);
    for id = 1:num_names
        if(metrics_1(id) > p_value_lim_log) || (metrics_2(id) > p_value_lim_log)
            names = vertcat(names, all_names(id));
            data_1 = vertcat(data_1, all_data_1(id, :));
            data_2 = vertcat(data_2, all_data_2(id, :));
        end
        
        if mod(id, 1000) == 0
            id = id
        end
    end
    
elseif strcmp(config.method, 'manova') && config.plot_method == 2

    p_value_lim = 1e-8;
    p_value_lim_log = -log10(p_value_lim);
    
    metrics_id = get_metrics_id(config);
    metrics_1 = all_data_1(:, metrics_id);
    metrics_2 = all_data_2(:, metrics_id);
    metrics_1 = process_metrics(metrics_1, config);
    metrics_2 = process_metrics(metrics_2, config);
    
    for id = 1:num_names
        if(metrics_1(id) < p_value_lim_log) || (metrics_2(id) < p_value_lim_log)
            names = vertcat(names, all_names(id));
            data_1 = vertcat(data_1, all_data_1(id, :));
            data_2 = vertcat(data_2, all_data_2(id, :));
        end
        
        if mod(id, 1000) == 0
            id = id
        end
    end
    
else
    
    names = all_names;
    data_1 = all_data_1;
    data_2 = all_data_2;
    
end



end