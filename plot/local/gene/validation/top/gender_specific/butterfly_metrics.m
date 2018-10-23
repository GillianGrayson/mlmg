clear all;

% ======== params ========
part = 0.05;
num_bins = 100;
diff_type = 2;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'gene_data';

config.chromosome_type = 'non_gender';

config.geo_type = 'islands_shores';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_variance';

config.disease = 'any';
config.gender = 'any';

config.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI') 
    config.up = 'D:/YandexDisk/Work/mlmg'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    config.up = 'D:/Aaron/Bio/mlmg'; 
else 
    config.up = 'E:/YandexDisk/Work/mlmg'; 
end 

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

if strcmp(getenv('computername'), 'MSI') 
    save_config.up = 'C:/Users/user/Google Drive/mlmg/figures'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    save_config.up = 'D:/Aaron/Bio/mlmg/figures'; 
else 
    save_config.up = 'C:/Users/user/Google Drive/mlmg/figures'; 
end
save_config.is_clustering = config.is_clustering;

% ======== processing ========
metrics_id_1 = get_metrics_id(config, 1);
metrics_id_2 = get_metrics_id(config, 2);

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
f_metrics_1 = f_top_data.data(:, metrics_id_1);
f_metrics_2 = f_top_data.data(:, metrics_id_2);
f_metrics_1 = process_metrics(f_metrics_1, config);
f_metrics_2 = process_metrics(f_metrics_2, config);

config.gender = 'M';
m_fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));
m_top_data = importdata(m_fn);
m_genes = m_top_data.textdata;
m_metrics_1 = m_top_data.data(:, metrics_id_1);
m_metrics_2 = m_top_data.data(:, metrics_id_2);
m_metrics_1 = process_metrics(m_metrics_1, config);
m_metrics_2 = process_metrics(m_metrics_2, config);

num_genes = size(f_genes, 1);

f_metrics_1_passed = f_metrics_1;
f_metrics_2_passed = f_metrics_2;
m_metrics_1_passed = zeros(num_genes, 1);
m_metrics_2_passed = zeros(num_genes, 1);
metrics_1_diff = zeros(num_genes, 1);
metrics_2_diff = zeros(num_genes, 1);
versus_metrics_diff = zeros(num_genes, 1);

for gene_id = 1:num_genes
    f_gene = f_genes(gene_id);
    m_id = find(m_genes==string(f_gene));
    m_metrics_1_passed(gene_id) = m_metrics_1(m_id);
    m_metrics_2_passed(gene_id) = m_metrics_2(m_id);
    metrics_1_diff(gene_id) = f_metrics_1_passed(gene_id) - m_metrics_1_passed(gene_id);
    metrics_2_diff(gene_id) = f_metrics_2_passed(gene_id) - m_metrics_2_passed(gene_id);
    versus_metrics_diff(gene_id) = process_diff(metrics_1_diff(gene_id), metrics_2_diff(gene_id), diff_type);
end

metrics_1_diff = normalize_metrics(metrics_1_diff, config);
metrics_2_diff = normalize_metrics(metrics_2_diff, config);

for gene_id = 1:num_genes
    versus_metrics_diff(gene_id) = process_diff(metrics_1_diff(gene_id), metrics_2_diff(gene_id), diff_type);
end

[versus_metrics_diff_srt, order] = sort(abs(versus_metrics_diff), 'descend');
genes_srt = f_genes;
metrics_1_diff_srt = zeros(num_genes, 1);
metrics_2_diff_srt = zeros(num_genes, 1);
for gene_id = 1:num_genes
    genes_srt(gene_id) = f_genes(order(gene_id));
    metrics_1_diff_srt(gene_id) = metrics_1_diff(order(gene_id));
    metrics_2_diff_srt(gene_id) = metrics_2_diff(order(gene_id));
end

num_rare = floor(part * num_genes);

common_genes = genes_srt(num_rare + 1:end);
common_metrics_1_diff = metrics_1_diff_srt(num_rare + 1:end);
common_metrics_2_diff = metrics_2_diff_srt(num_rare + 1:end);

rare_genes = genes_srt(1:num_rare);
rare_metrics_1_diff = metrics_1_diff_srt(1:num_rare);
rare_metrics_2_diff = metrics_2_diff_srt(1:num_rare);

f1 = figure;
hold all;
h = plot([min(metrics_1_diff) max(metrics_1_diff)], [0 0], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(metrics_2_diff) max(metrics_2_diff)], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all
h = plot(common_metrics_1_diff, common_metrics_2_diff, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 30);
xlabel('$\Delta$ metrics 1', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$ metrics 2', 'Interpreter', 'latex');

for gene_id = 1:size(rare_genes, 1)
    hold all;
    h = plot(rare_metrics_1_diff(gene_id), rare_metrics_2_diff(gene_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, string(rare_genes(gene_id)))
end

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f1, sprintf('%s/butterfly_metrics_scatter_%s.fig', save_path, suffix))
saveas(f1, sprintf('%s/butterfly_metrics_scatter_%s.png', save_path, suffix))

metrics_1_diff_begin = min(metrics_1_diff_srt);
metrics_1_diff_end = max(metrics_1_diff_srt);
metrics_1_diff_step = (metrics_1_diff_end - metrics_1_diff_begin) / num_bins;
metrics_1_diff_bins = linspace(metrics_1_diff_begin +  0.5 * metrics_1_diff_step, metrics_1_diff_end - 0.5 * metrics_1_diff_step, num_bins);

metrics_2_diff_begin = min(metrics_2_diff_srt);
metrics_2_diff_end = max(metrics_2_diff_srt);
metrics_2_diff_step = (metrics_2_diff_end - metrics_2_diff_begin) / num_bins;
metrics_2_diff_bins = linspace(metrics_2_diff_begin +  0.5 * metrics_2_diff_step, metrics_2_diff_end - 0.5 * metrics_2_diff_step, num_bins);

metrics_pdf = zeros(num_bins);
for gene_id = 1:num_genes
    id = floor((metrics_1_diff_srt(gene_id) - metrics_1_diff_begin) * num_bins / (metrics_1_diff_end - metrics_1_diff_begin + 1e-8)) + 1;
    y_id = floor((metrics_2_diff_srt(gene_id) - metrics_2_diff_begin) * num_bins / (metrics_2_diff_end - metrics_2_diff_begin + 1e-8)) + 1;
    metrics_pdf(id, y_id) = metrics_pdf(id, y_id) + 1;
end

metrics_pdf = metrics_pdf / (sum(sum(metrics_pdf)) * metrics_1_diff_step * metrics_2_diff_step);
norm = sum(sum(metrics_pdf)) * metrics_1_diff_step * metrics_2_diff_step

f2 = figure;
h = imagesc(metrics_1_diff_bins, metrics_2_diff_bins, metrics_pdf');
set(gca, 'FontSize', 30);
xlabel('$\Delta$ metrics 1', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$ metrics 2', 'Interpreter', 'latex');
colormap hot;
cb = colorbar;
set(gca, 'FontSize', 30);
title(cb, 'PDF');
set(gca,'YDir','normal');

hold all;
h = plot([min(metrics_1_diff_srt) max(metrics_1_diff_srt)], [0 0], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(metrics_2_diff_srt) max(metrics_2_diff_srt)], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

for gene_id = 1:size(rare_genes, 1)
    hold all;
    h = plot(rare_metrics_1_diff(gene_id), rare_metrics_2_diff(gene_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, string(rare_genes(gene_id)))
end

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f2, sprintf('%s/butterfly_metrics_pdf_%s.fig', save_path, suffix))
saveas(f2, sprintf('%s/butterfly_metrics_pdf_%s.png', save_path, suffix))

abs_diff = abs(versus_metrics_diff_srt);
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

savefig(f3, sprintf('%s/butterfly_metrics_delta_%s.fig', save_path, suffix))
saveas(f3, sprintf('%s/butterfly_metrics_delta_%s.png', save_path, suffix))