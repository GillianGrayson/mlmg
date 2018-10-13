clear all;

% ======== params ========
part = 0.05;
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
metrics_label = get_metrics_label(config);

save_path = sprintf('%s/%s', ...
    save_config.up, ...
    get_result_path(save_config));
mkdir(save_path);
suffix = sprintf('method(%s)_part(%0.4f)', ...
    config.method, ...
    part);

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
metrics_diff = zeros(num_genes, 1);

for gene_id = 1:num_genes
    f_gene = f_genes(gene_id);
    m_id = find(m_genes==string(f_gene));
    m_metrics_passed(gene_id) = m_metrics(m_id);
    metrics_diff(gene_id) = f_metrics_passed(gene_id) - m_metrics_passed(gene_id);
end

[metrics_diff_srt, order] = sort(abs(metrics_diff), 'descend');
genes_srt = f_genes;
f_ids_srt = zeros(num_genes, 1);
m_ids_srt = zeros(num_genes, 1);
f_metrics_srt = zeros(num_genes, 1);
m_metrics_srt = zeros(num_genes, 1);
for gene_id = 1:num_genes
    genes_srt(gene_id) = f_genes(order(gene_id));
    f_metrics_srt(gene_id) = f_metrics_passed(order(gene_id));
    m_metrics_srt(gene_id) = m_metrics_passed(order(gene_id));
end

num_rare = floor(part * num_genes);

common_genes = genes_srt(num_rare + 1:end);
common_f_metrics = f_metrics_srt(num_rare + 1:end);
common_m_metrics = m_metrics_srt(num_rare + 1:end);
common_diff = metrics_diff_srt(num_rare + 1:end);

rare_genes = genes_srt(1:num_rare);
rare_f_metrics = f_metrics_srt(1:num_rare);
rare_m_metrics = m_metrics_srt(1:num_rare);
rare_diff = metrics_diff_srt(1:num_rare);

f1 = figure;
hold all;
h = plot([min(f_metrics_srt) max(f_metrics_srt)], [0 0], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(m_metrics_srt) max(m_metrics_srt)], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all
h = plot(common_f_metrics, common_m_metrics, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 30);
xlabel(sprintf('%s F', metrics_label), 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(sprintf('%s M', metrics_label), 'Interpreter', 'latex');

for gene_id = 1:size(rare_genes, 1)
    hold all;
    h = plot(rare_f_metrics(gene_id), rare_m_metrics(gene_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, string(rare_genes(gene_id)))
end

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f1, sprintf('%s/butterfly_scatter_%s.fig', save_path, suffix))
saveas(f1, sprintf('%s/butterfly_scatter_%s.png', save_path, suffix))

f_metrics_begin = min(f_metrics_srt);
f_metrics_end = max(f_metrics_srt);
f_metrics_step = (f_metrics_end - f_metrics_begin) / num_bins;
f_metrics_bins = linspace(f_metrics_begin +  0.5 * f_metrics_step, f_metrics_end - 0.5 * f_metrics_step, num_bins);

m_metrics_begin = min(m_metrics_srt);
m_metrics_end = max(m_metrics_srt);
m_metrics_step = (m_metrics_end - m_metrics_begin) / num_bins;
m_metrics_bins = linspace(m_metrics_begin +  0.5 * m_metrics_step, m_metrics_end - 0.5 * m_metrics_step, num_bins);

metrics_pdf = zeros(num_bins);
for gene_id = 1:num_genes
    id = floor((f_metrics_srt(gene_id) - f_metrics_begin) * num_bins / (f_metrics_end - f_metrics_begin + 1e-8)) + 1;
    y_id = floor((m_metrics_srt(gene_id) - m_metrics_begin) * num_bins / (m_metrics_end - m_metrics_begin + 1e-8)) + 1;
    metrics_pdf(id, y_id) = metrics_pdf(id, y_id) + 1;
end

metrics_pdf = metrics_pdf / (sum(sum(metrics_pdf)) * f_metrics_step * m_metrics_step);
norm = sum(sum(metrics_pdf)) * f_metrics_step * m_metrics_step

f2 = figure;
h = imagesc(f_metrics_bins, m_metrics_bins, metrics_pdf');
set(gca, 'FontSize', 30);
xlabel(sprintf('%s F', metrics_label), 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(sprintf('%s M', metrics_label), 'Interpreter', 'latex');
colormap hot;
cb = colorbar;
set(gca, 'FontSize', 30);
title(cb, 'PDF');
set(gca,'YDir','normal');

hold all;
h = plot([min(f_metrics_srt) max(f_metrics_srt)], [0 0], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(m_metrics_srt) max(m_metrics_srt)], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

for gene_id = 1:size(rare_genes, 1)
    hold all;
    h = plot(rare_f_metrics(gene_id), rare_m_metrics(gene_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, string(rare_genes(gene_id)))
end

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f2, sprintf('%s/butterfly_pdf_%s.fig', save_path, suffix))
saveas(f2, sprintf('%s/butterfly_pdf_%s.png', save_path, suffix))

abs_diff = abs(metrics_diff_srt);
diff_begin = min(abs_diff);
diff_end = max(abs_diff);
diff_step = (diff_end - diff_begin) / num_bins;
diff_bins = linspace(diff_begin +  0.5 * diff_step, diff_end - 0.5 * diff_step, num_bins);

diff_pdf = zeros(num_bins, 1);
for gene_id = 1:num_genes
    id = floor((abs_diff(gene_id) - diff_begin) * num_bins / (diff_end - diff_begin + 1e-8)) + 1;
    diff_pdf(id) = diff_pdf(id) + 1;
end
diff_pdf = diff_pdf / (sum(diff_pdf) * diff_step);
norm = sum(diff_pdf) * diff_step

f3 = figure;
hold all
h = plot(diff_bins, diff_pdf, 'LineWidth', 2, 'Color', 'r');
set(gca, 'FontSize', 30);
xlabel('$|\Delta|$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');

limit_id = floor((abs_diff(num_rare) - diff_begin) * num_bins / (diff_end - diff_begin + 1e-8)) + 1;
hold all;
h = plot(diff_bins(limit_id:end), diff_pdf(limit_id:end), 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(end)], [0 0], 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(limit_id)], [0 diff_pdf(limit_id)], 'LineWidth', 4, 'Color', 'k');

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f3, sprintf('%s/delta_pdf_%s.fig', save_path, suffix))
saveas(f3, sprintf('%s/delta_pdf_%s.png', save_path, suffix))