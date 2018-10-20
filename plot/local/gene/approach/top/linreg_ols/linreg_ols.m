clear all;

% ======== params ========
gene = 'TLE1';

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'gene_data';

config.chromosome_type = 'non_gender';

config.geo_type = 'islands_shores';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_ols';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI')
    config.up = 'D:/YandexDisk/Work/mlmg';
else
    config.up = 'E:/YandexDisk/Work/mlmg';
end

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    plot_linreg_ols(config, gene)
    config.gender = 'M';
    plot_linreg_ols(config, gene)
    config.gender = 'versus';
else
    plot_linreg_ols(config, gene)
end

suffix = sprintf('gene(%s)_gender(%s)', gene, config.gender);


up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/linreg_ols_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/linreg_ols_%s.png', save_path, suffix))


function plot_linreg_ols(config, gene)

fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));

top_data = importdata(fn);

genes = top_data.textdata;
intercepts = top_data.data(:, 2);
slopes = top_data.data(:, 3);
intercepts_std = top_data.data(:, 4);
slopes_std = top_data.data(:, 5);

indexes = get_attributes_indexes(config);
ages = get_ages(config);

ages_passed = zeros(size(indexes, 1), 1);
for id = 1:size(indexes, 1)
    index = indexes(id);
    ages_passed(id) = ages(index);
end

fn = sprintf('%s/data/%s/gene_data.txt', ...
    config.up, ...
    get_gene_data_path(config));

data = importdata(fn);
genes_names = data.textdata;
genes_data = data.data;

gene_name = string(gene);
idx = find(genes_names==gene_name);
gene_data = genes_data(idx, :)';

gene_data_passed = size(indexes, 1);
for id = 1:size(indexes, 1)
    gene_data_passed(id) = gene_data(indexes(id));
end

gene_id = find(genes==gene_name);

slope = slopes(gene_id);
intercept = intercepts(gene_id);
intercept_std = intercepts_std(gene_id);
slope_std = slopes_std(gene_id);
intercept_down = intercept - 3 * intercept_std;
slope_down = slope - 3 * slope_std;
intercept_up = intercept + 3 * intercept_std;
slope_up = slope + 3 * slope_std;
x_lin = [min(ages), max(ages)];
y_lin = [slope * x_lin(1) + intercept, slope * x_lin(2) + intercept];
y_down = [slope_down * x_lin(1) + intercept_down, slope_down * x_lin(2) + intercept_down];
y_up = [slope_up * x_lin(1) + intercept_up, slope_up * x_lin(2) + intercept_up];

hold all;
h = plot(ages_passed, gene_data_passed, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
color = get(h, 'Color');

hold all;
h = plot(x_lin, y_lin, '-', 'LineWidth', 3);
legend(h, sprintf('%s: %s', gene_name, config.gender));
set(h, 'Color', color)
set(gca, 'FontSize', 30);
xlabel('age', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\beta$', 'Interpreter', 'latex');
xlim([min(ages) - (max(ages) - min(ages)) * 0.1, max(ages) + (max(ages) - min(ages)) * 0.1])

hold all;
h = plot(x_lin, y_down, '-.', 'LineWidth', 1, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot(x_lin, y_up, '-.', 'LineWidth', 1, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

box on;
end






