clear all;

age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

base = 'GSE87571';
data_type = 'gene_data';

chromosome_type = 'non_gender';

geo_type = 'islands_shores';
gene_data_type = 'mean';

info_type = 'result';

disease = 'any';
scenario = 'approach';
approach = 'top';
method = 'linreg';

rare_part = 0.005;

metrics_id = 1;
if strcmp(method, 'linreg')
    metrics_id = 3;
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
f_genes = f_top_data.textdata;
f_metrics = f_top_data.data(:, metrics_id);

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
m_genes = m_top_data.textdata;
m_metrics = m_top_data.data(:, metrics_id);

num_genes = size(f_genes, 1);

f_ids = linspace(1, num_genes, num_genes);
m_ids = zeros(num_genes, 1);

f_metrics_passed = f_metrics;
m_metrics_passed = zeros(num_genes, 1);

metrics_diff = zeros(num_genes, 1);

for gene_id = 1:num_genes
    f_gene = f_genes(gene_id);
    m_id = find(m_genes==string(f_gene));
    m_ids(gene_id) = m_id;
    m_metrics_passed(gene_id) = m_metrics(m_id);
    metrics_diff(gene_id) = f_metrics_passed(gene_id) - m_metrics_passed(gene_id);
end

[metrix_diff_srt, order] = sort(abs(metrics_diff), 'descend');
genes_srt = f_genes;
f_ids_srt = zeros(num_genes, 1);
m_ids_srt = zeros(num_genes, 1);
f_metrics_srt = zeros(num_genes, 1);
m_metrics_srt = zeros(num_genes, 1);
for gene_id = 1:num_genes
    genes_srt(gene_id) = f_genes(order(gene_id));
    f_ids_srt(gene_id) = f_ids(order(gene_id));
    m_ids_srt(gene_id) = m_ids(order(gene_id));
    f_metrics_srt(gene_id) = f_metrics_passed(order(gene_id));
    m_metrics_srt(gene_id) = m_metrics_passed(order(gene_id));
end

num_rare = floor(rare_part * num_genes);

common_genes = genes_srt(num_rare + 1:end);
common_f_metrics = f_metrics_srt(num_rare + 1:end);
common_m_metrics = m_metrics_srt(num_rare + 1:end);
common_diff = metrix_diff_srt(num_rare + 1:end);

rare_genes = genes_srt(1:num_rare);
rare_f_metrics = f_metrics_srt(1:num_rare);
rare_m_metrics = m_metrics_srt(1:num_rare);
rare_diff = metrix_diff_srt(1:num_rare);

figure
h = plot(common_f_metrics, common_m_metrics, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 30);
xlabel('metrics F', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('metrics M', 'Interpreter', 'latex');

for  gene_id = 1:size(rare_genes, 1)
    hold all;
    h = plot(rare_f_metrics(gene_id), rare_m_metrics(gene_id), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    legend(h, string(rare_genes(gene_id)))
end

hold all;
h = plot([min(f_metrics_srt) max(f_metrics_srt)], [0 0], 'LineWidth', 3, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(m_metrics_srt) max(m_metrics_srt)], 'LineWidth', 3, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

propertyeditor('on')
box on;
