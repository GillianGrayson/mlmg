clear all;

% ======== params ========
gene = 'FIGN';

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'gene_data';

config.cross_reactive = 'cross_reactive_excluded';
config.snp = 'snp_excluded_weak';

config.chromosome_type = 'non_gender';
config.geo_type = 'any';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'inside_gene';
config.method = 'moment';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI') 
    config.up = 'D:/YandexDisk/Work/mlmg'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    config.up = 'C:/Users/User/YandexDisk/mlmg'; 
else 
    config.up = 'E:/YandexDisk/Work/mlmg'; 
end 

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    plot_moment(config, gene)
    config.gender = 'M';
    plot_moment(config, gene)
    config.gender = 'versus';
else
    plot_moment(config, gene)
end

suffix = sprintf('gene(%s)_gender(%s)', gene, config.gender);

legend('-DynamicLegend');

if strcmp(getenv('computername'), 'MSI') 
    up_save = 'C:/Users/user/Google Drive/mlmg/figures'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    up_save = 'D:/Aaron/Bio/mlmg/figures'; 
else 
    up_save = 'C:/Users/user/Google Drive/mlmg/figures'; 
end 

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

savefig(f, sprintf('%s/moment_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/moment_%s.png', save_path, suffix))


function plot_moment(config, gene)
print_rate = 1000;
fn = sprintf('%s/data/%s/inside_gene.txt', ...
    config.up, ...
    get_result_path(config));

inside_gene_data = importdata(fn, ' ');
genes = inside_gene_data.textdata(:, 1);
cpgs = inside_gene_data.textdata(:, 2);
means = inside_gene_data.data(:, 1);
stds = inside_gene_data.data(:, 2);

ids = find(string(genes) == string(gene));
num_targets = size(ids, 1);
target_cpgs = {};
target_x = zeros(num_targets, 1);
target_y = zeros(num_targets, 1);
target_errors = zeros(num_targets, 1);
for i = 1:num_targets
    target_cpgs{i} = string(cpgs(ids(i)));
    target_x(i) = i;
    target_y(i) = means(ids(i));
    target_errors(i) = stds(ids(i));
end

hold all;

h = errorbar(target_x, target_y, target_errors, 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'LineWidth', 2, 'CapSize', 15);
legend(h, sprintf('%s : %s', gene, config.gender));
xticks(target_x);
xticklabels(target_cpgs);
xtickangle(45);
xlim([0, num_targets + 1]);
set(gca, 'FontSize', 30);
xL = xlabel('', 'Interpreter', 'latex');
yL = ylabel('$\beta$', 'Interpreter', 'latex');

ax = ancestor(h, 'axes');
xrule = ax.XAxis;
xrule.FontSize = 15;

propertyeditor('on')
box on;
end

