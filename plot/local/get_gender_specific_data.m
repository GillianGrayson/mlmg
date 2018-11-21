function [names, data_1, data_2] = get_gender_specific_data(config)

suffix = '';
if isfield(config, 'suffix')
    suffix = config.suffix;
end

config.gender = 'F';
fn = sprintf('%s/data/%s/top%s.txt', ...
    config.up, ...
    get_result_path(config), ...
    suffix);
raw_data_1 = importdata(fn, ' ');

config.gender = 'M';
fn = sprintf('%s/data/%s/top%s.txt', ...
    config.up, ...
    get_result_path(config), ...
    suffix);
raw_data_2 = importdata(fn, ' ');

names_1 = raw_data_1.textdata;
d_1_tmp = raw_data_1.data;

map_1 = containers.Map();
for id = 1:size(names_1, 1)
    map_1(string(names_1(id))) = d_1_tmp(id, :);
end

names_2 = raw_data_2.textdata;
d_2_tmp = raw_data_2.data;

map_2 = containers.Map();
for id = 1:size(names_1, 1)
    map_2(string(names_2(id))) = d_2_tmp(id, :);
end

all_names = names_1;
num_names = size(all_names, 1);

all_data_1 = d_1_tmp;
all_data_2 = zeros(num_names, 1);
for id_1 = 1:num_names
    name = string(all_names(id_1));
    data_2_curr = map_2(name);
    for m_id = 1:size(all_data_1, 2)
        all_data_2(id_1, m_id) = data_2_curr(m_id);
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
    
elseif strcmp(config.method, 'manova') && config.plot_method == 3
    
    p_value_lim = 1e-8;
    p_value_lim_log = -log10(p_value_lim);
    
    metrics_id = get_metrics_id(config);
    metrics_1 = all_data_1(:, metrics_id);
    metrics_2 = all_data_2(:, metrics_id);
    metrics_1 = process_metrics(metrics_1, config);
    metrics_2 = process_metrics(metrics_2, config);
    
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
    
elseif strcmp(config.method, 'manova') && config.plot_method == 4
    metrics_id = get_metrics_id(config);
    metrics_1 = all_data_1(:, metrics_id);
    metrics_2 = all_data_2(:, metrics_id);
    metrics_1 = process_metrics(metrics_1, config);
    metrics_2 = process_metrics(metrics_2, config);
    
    config_moment.data_base = config.data_base;
    config_moment.data_type = config.data_type;
    
    config_moment.chromosome_type = config.chromosome_type;
    
    config_moment.class_type = config.class_type;
    
    config_moment.info_type = config.info_type;
    
    config_moment.scenario = config.scenario;
    config_moment.approach = 'inside_bop';
    config_moment.method = 'moment';
    
    config_moment.disease = 'any';
    config_moment.gender = 'F';
    
    config_moment.is_clustering = 0;
    
    if strcmp(getenv('computername'), 'MSI')
        config_moment.up = 'D:/YandexDisk/Work/mlmg';
    elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS')
        config_moment.up = 'D:/Aaron/Bio/mlmg';
    else
        config_moment.up = 'E:/YandexDisk/Work/mlmg';
    end
    
    fn = sprintf('%s/data/%s/inside_bop.txt', ...
        config_moment.up, ...
        get_result_path(config_moment));
    
    f_inside_bop_data = importdata(fn, ' ');
    f_bops = f_inside_bop_data.textdata(:, 1);
    f_means = f_inside_bop_data.data(:, 1);
    f_stds = f_inside_bop_data.data(:, 2);
    
    config_moment.gender = 'M';
    
    fn = sprintf('%s/data/%s/inside_bop.txt', ...
        config_moment.up, ...
        get_result_path(config_moment));
    
    m_inside_bop_data = importdata(fn, ' ');
    m_bops = m_inside_bop_data.textdata(:, 1);
    m_means = m_inside_bop_data.data(:, 1);
    m_stds = m_inside_bop_data.data(:, 2);
    
    for id = 1:size(metrics_1, 1)
        bop_name = all_names{id};
        ids = find(string(f_bops) == string(bop_name));
        num_targets = size(ids, 1);
        for i = 1:num_targets
            f_interval = [f_means(ids(i))-f_stds(ids(i)), f_means(ids(i))+f_stds(ids(i))];
            m_interval = [m_means(ids(i))-m_stds(ids(i)), m_means(ids(i))+m_stds(ids(i))];
            if((f_interval(1)<min(m_interval)) && (min(m_interval)<f_interval(2))) || ((f_interval(1)<max(m_interval)) && (max(m_interval)<f_interval(2))) || ((m_interval(1)<min(f_interval)) && (min(f_interval)<m_interval(2))) || ((m_interval(1)<max(f_interval)) && (max(f_interval)<m_interval(2)))
                continue
            else
                names = vertcat(names, all_names(id));
                data_1 = vertcat(data_1, all_data_1(id, :));
                data_2 = vertcat(data_2, all_data_2(id, :));
            end
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