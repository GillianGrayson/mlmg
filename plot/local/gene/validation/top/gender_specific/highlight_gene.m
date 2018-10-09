clear all;

age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

base = 'data_base_versus';
data_type = 'gene_data';

chromosome_type = 'non_gender';

geo_type = 'islands_shores';
gene_data_type = 'mean';

info_type = 'result';

disease = 'any';
scenario = 'validation';
approach = 'top';
method = 'gender_specific';

target_data_bases = 'GSE40279_GSE87571';
target_method = 'linreg';
target_part = 0.05;

up = '../../../../../..';

fn = sprintf('%s/data/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/intersection_genes_data_bases(%s)_method(%s)_part(%0.2f).txt', ...
    up, ...
    base, ...
    data_type, ...
    chromosome_type, ...
    geo_type, ...
    gene_data_type, ...
    info_type, ...
    disease, ...
    'F', ...
    scenario, ...
    approach, ...
    method, ...
    target_data_bases, ...
    target_method, ...
    target_part);
target_genes = importdata(fn);

base = 'GSE40279';
data_type = 'gene_data';

chromosome_type = 'non_gender';

geo_type = 'islands_shores';
gene_data_type = 'mean';

info_type = 'result';

disease = 'any';
scenario = 'approach';
approach = 'top';
method = 'linreg';

is_clustering = 0;

rare_part = 0.005;

num_bins = 100;

metrics_id = 1;
if strcmp(method, 'linreg')
    if is_clustering == 1
        metrics_id = 3;
    else
        metrics_id = 1;
    end
end

up = '../../../../../..';

f_fn = sprintf('%s/data/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/top.txt', ...
    up, ...
    base, ...
    data_type, ...
    chromosome_type, ...
    geo_type, ...
    gene_data_type, ...
    info_type, ...
    disease, ...
    'F', ...
    scenario, ...
    approach, ...
    method);
f_top_data = importdata(f_fn);
f_genes = f_top_data.textdata(:, 1);
f_metrics = f_top_data.data(:, metrics_id);

if strcmp(method, 'linreg')
    nothing_to_do = 0;
elseif strcmp(method, 'manova')
    f_metrics = -log10(f_metrics);
end

m_fn = sprintf('%s/data/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/top.txt', ...
    up, ...
    base, ...
    data_type, ...
    chromosome_type, ...
    geo_type, ...
    gene_data_type, ...
    info_type, ...
    disease, ...
    'M', ...
    scenario, ...
    approach, ...
    method);
m_top_data = importdata(m_fn);
m_genes = m_top_data.textdata(:, 1);
m_metrics = m_top_data.data(:, metrics_id);

if strcmp(method, 'linreg')
    nothing_to_do = 0;
elseif strcmp(method, 'manova')
    m_metrics = -log10(m_metrics);
end

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


