function [names, metrics_1, metrics_2] = process_gender_specific_metrics(config)
names = [];
metrics_1 = [];
metrics_2 = [];
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
    
    config.gender = 'F';
    f_fn = sprintf('%s/data/%s/top.txt', ...
        config.up, ...
        get_result_path(config));
    f_top_data = importdata(f_fn);
    f_genes = f_top_data.textdata;
    f_metrics = f_top_data.data(:, metrics_id);
    f_metrics = process_metrics(f_metrics, config);
    
    config.gender = 'M';
    m_fn = sprintf('%s/data/%s/top.txt', ...
        config.up, ...
        get_result_path(config));
    m_top_data = importdata(m_fn);
    m_genes = m_top_data.textdata;
    m_metrics = m_top_data.data(:, metrics_id);
    m_metrics = process_metrics(m_metrics, config);
    
    num_genes = size(f_genes, 1);
    
    f_metrics_passed = f_metrics;
    m_metrics_passed = zeros(num_genes, 1);
    for gene_id = 1:num_genes
        f_gene = f_genes(gene_id);
        m_id = find(m_genes==string(f_gene));
        m_metrics_passed(gene_id) = m_metrics(m_id);
    end
    
    [f_metrics_passed, m_metrics_passed] = process_metrics_plane(f_metrics_passed, m_metrics_passed, config);
    num_bops = size(f_metrics_passed, 1);
    
    
    
    
    
    
    
    
    
    
    
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