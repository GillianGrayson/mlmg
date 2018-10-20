clear all;

% ======== params ========
gene = 'NCRNA00181';

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
    plot_linreg_variance(config, gene)
    config.gender = 'M';
    plot_linreg_variance(config, gene)
    config.gender = 'versus';
else
    plot_linreg_variance(config, gene)
end

suffix = sprintf('gene(%s)_gender(%s)', gene, config.gender);


up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/linreg_variance_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/linreg_variance_%s.png', save_path, suffix))


function plot_linreg_variance(config, gene)

fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));

top_data = importdata(fn);

genes = top_data.textdata;
slopes = top_data.data(:, 3);
intercepts = top_data.data(:, 4);
slopes_variance = top_data.data(:, 7);
intercepts_variance = top_data.data(:, 8);

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
slope_variance = slopes_variance(gene_id);
intercept_variance = intercepts_variance(gene_id);
x_lin = [min(ages), max(ages)];
y_lin = [slope * x_lin(1) + intercept, slope * x_lin(2) + intercept];

errors = zeros(size(ages_passed, 1), 1);
for i = 1:size(ages_passed, 1)
    errors(i) = abs(gene_data_passed(i) - (ages_passed(i) * slope + intercept));
end
errors_lin_pos = [slope_variance * x_lin(1) + intercept_variance, slope_variance * x_lin(2) + intercept_variance];
errors_lin_neg = [-slope_variance * x_lin(1) - intercept_variance, -slope_variance * x_lin(2) - intercept_variance];

pos_p1 = [x_lin(1), errors_lin_pos(1) + intercept]';
pos_p2 = [x_lin(2), errors_lin_pos(2) + intercept]';
angle_pos = atan(slope);
rotate_mtx_pos = [cos(angle_pos) -sin(angle_pos); sin(angle_pos) cos(angle_pos)];

neg_p1 = [x_lin(1), errors_lin_neg(1) + intercept]';
neg_p2 = [x_lin(2), errors_lin_neg(2) + intercept]';

var_pos_p1 = rotate_mtx_pos * pos_p1;
var_pos_p2 = rotate_mtx_pos * pos_p2;

var_neg_p1 = rotate_mtx_pos * neg_p1;
var_neg_p2 = rotate_mtx_pos * neg_p2;

hold all;
h = plot(ages_passed, gene_data_passed, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
color = get(h, 'Color');

hold all;
h = plot(x_lin, y_lin, '-', 'LineWidth', 3);
legend(h, sprintf('%s: %s', gene_name, config.gender));
set(h, 'Color', color)

hold all;
h = plot([var_pos_p1(1),  var_pos_p2(1)], [var_pos_p1(2),  var_pos_p2(2)], '-', 'LineWidth', 1, 'Color', color);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([var_neg_p1(1),  var_neg_p2(1)], [var_neg_p1(2),  var_neg_p2(2)], '-', 'LineWidth', 1, 'Color', color);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'FontSize', 30);
xlabel('age', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\beta$', 'Interpreter', 'latex');
xlim([min(ages) - (max(ages) - min(ages)) * 0.1, max(ages) + (max(ages) - min(ages)) * 0.1])
box on;
end
