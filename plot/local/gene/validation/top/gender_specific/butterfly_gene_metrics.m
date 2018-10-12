clear all;

% ======== params ========
num_bins = 100;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'gene_data';

config.chromosome_type = 'non_gender';

config.geo_type = 'islands_shores';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg';

config.disease = 'any';
config.gender = 'any';

config.up = '../../../../../..';
config.is_clustering = 0;

% ======== save_config ========
save_config.data_base = config.data_base;
save_config.data_type = config.data_type;

save_config.chromosome_type = config.chromosome_type;

save_config.geo_type = config.geo_type;
save_config.gene_data_type = config.gene_data_type;

save_config.info_type = 'result';

save_config.scenario = 'validation';
save_config.approach = 'top';
save_config.method = 'gender_specific';

save_config.disease = config.disease;
save_config.gender = 'versus';

save_config.up = 'C:/Users/user/Google Drive/mlmg/figures';
save_config.is_clustering = config.is_clustering;

% ======== processing ========
metrics_id = get_metrics_id(config);
add_metrics_id = get_add_metrics_id(config);

save_path = sprintf('%s/%s', ...
    save_config.up, ...
    get_result_path(save_config));
mkdir(save_path);
suffix = sprintf('method(%s)', ...
    config.method);

config.gender = 'F';
f_fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));
f_top_data = importdata(f_fn);
f_genes = f_top_data.textdata;
f_metrics = f_top_data.data(:, metrics_id);
f_add_metrics = f_top_data.data(:, add_metrics_id);
f_metrics = process_metrics(f_metrics, config);
f_add_metrics = process_metrics(f_add_metrics, config);

config.gender = 'M';
m_fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));
m_top_data = importdata(m_fn);
m_genes = m_top_data.textdata;
m_metrics = m_top_data.data(:, metrics_id);
m_add_metrics = m_top_data.data(:, add_metrics_id);
m_metrics = process_metrics(m_metrics, config);
m_add_metrics = process_metrics(m_add_metrics, config);

num_genes = size(f_genes, 1);

f_metrics_passed = f_metrics;
f_add_metrics_passed = f_add_metrics;
m_metrics_passed = zeros(num_genes, 1);
m_add_metrics_passed = zeros(num_genes, 1);
metrics_diff = zeros(num_genes, 1);
add_metrics_diff = zeros(num_genes, 1);

for gene_id = 1:num_genes
    f_gene = f_genes(gene_id);
    m_id = find(m_genes==string(f_gene));
    m_metrics_passed(gene_id) = m_metrics(m_id);
    m_add_metrics_passed(gene_id) = m_add_metrics(m_id);
    metrics_diff(gene_id) = f_metrics_passed(gene_id) - m_metrics_passed(gene_id);
    add_metrics_diff(gene_id) = f_add_metrics_passed(gene_id) - m_add_metrics_passed(gene_id);
end

f1 = figure;
hold all;
h = plot([min(metrics_diff) max(metrics_diff)], [0 0], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(add_metrics_diff) max(add_metrics_diff)], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all
h = plot(metrics_diff, add_metrics_diff, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 30);
xlabel('$\Delta$ metrics 1', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$ metrics 2', 'Interpreter', 'latex');

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f1, sprintf('%s/butterfly_metrics_scatter_%s.fig', save_path, suffix))
saveas(f1, sprintf('%s/butterfly_metrics_scatter_%s.png', save_path, suffix))