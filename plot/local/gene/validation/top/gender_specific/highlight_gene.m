clear all;

% ======== params ========
target_data_bases = 'GSE40279_GSE87571';
target_method = 'linreg';
target_part = 0.05;

% ======== config ========
config.data_base = 'data_base_versus';
config.data_type = 'gene_data';

config.chromosome_type = 'non_gender';

config.geo_type = 'islands_shores';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'validation';
config.approach = 'top';
config.method = 'gender_specific';

config.disease = 'any';
config.gender = 'any';

config.up = '../../../../../..';
config.is_clustering = 0;

% ======== target_config ========
target_config.data_base = 'GSE40279';
target_config.data_type = config.data_type;

target_config.chromosome_type = config.chromosome_type;

target_config.geo_type = config.geo_type;
target_config.gene_data_type = config.gene_data_type;

target_config.info_type = config.info_type;

target_config.scenario = 'approach';
target_config.approach = 'top';
target_config.method = 'linreg';

target_config.disease = config.disease;
target_config.gender = config.gender;

target_config.up = config.up;
target_config.is_clustering = config.is_clustering;

% ======== processing ========
suffix = sprintf('data_bases(%s)_method(%s)_part(%0.2f)', ...
    target_data_bases, ...
    target_method, ...
    target_part);

fn = sprintf('%s/data/%s/intersection_genes_data_%s.txt', ...
    config.up, ...
    get_result_path(config), ...
    suffix);
target_genes = importdata(fn);

metrics_id = get_metrics_id(config);

target_config.gender = 'F';
f_fn = sprintf('%s/data/%s/top.txt', ...
    target_config.up, ...
    get_result_path(target_config));
f_top_data = importdata(f_fn);
f_genes = f_top_data.textdata;
f_metrics = f_top_data.data(:, metrics_id);
f_metrics = process_metrics(f_metrics, target_config);

target_config.gender = 'M';
m_fn = sprintf('%s/data/%s/top.txt', ...
    target_config.up, ...
    get_result_path(target_config));
m_top_data = importdata(m_fn);
m_genes = m_top_data.textdata;
m_metrics = m_top_data.data(:, metrics_id);
m_metrics = process_metrics(m_metrics, target_config);

num_genes = size(target_genes, 1);

target_f_metrics = zeros(num_genes, 1);
target_m_metrics = zeros(num_genes, 1);
target_diff = zeros(num_genes, 1);

for gene_id = 1:num_genes
    gene = target_genes(gene_id);
    f_id = find(f_genes==string(gene));
    m_id = find(m_genes==string(gene));
    target_f_metrics(gene_id) = f_metrics(f_id);
    target_m_metrics(gene_id) = m_metrics(m_id);
    target_diff(gene_id) = target_f_metrics(gene_id) - target_m_metrics(gene_id);
end

for gene_id = 1:num_genes
    hold all;
    h = plot(target_f_metrics(gene_id), target_m_metrics(gene_id), 'Marker', 'square', 'MarkerSize', 15, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, sprintf('%s', string(target_genes(gene_id))))
end

propertyeditor('on')
box on;


