clear all;

% ======== params ========
bop = 'chr20:44657463-44659243*Island';

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'bop_data';

config.chromosome_type = 'non_gender';

config.class_type = 'ClassAB';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'inside_bop';
config.method = 'moment';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;
config.up = '../../../../../..';

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    plot_moment(config, bop)
    config.gender = 'M';
    plot_moment(config, bop)
    config.gender = 'versus';
else
    plot_moment(config, bop)
end

bop_name = strrep(bop, ':', '_');
bop_name = strrep(bop_name, '*', '_');
bop_name = strrep(bop_name, '-', '_');
suffix = sprintf('bop(%s)_gender(%s)', bop_name, config.gender);

legend('-DynamicLegend');

up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

savefig(f, sprintf('%s/moment_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/moment_%s.png', save_path, suffix))


function plot_moment(config, bop)
print_rate = 1000;
fn = sprintf('%s/data/%s/inside_bop.txt', ...
    config.up, ...
    get_result_path(config));

inside_bop_data = importdata(fn, ' ');
bops = inside_bop_data.textdata(:, 1);
genes = inside_bop_data.textdata(:, 2);
cpgs = inside_bop_data.textdata(:, 3);
means = inside_bop_data.data(:, 1);
stds = inside_bop_data.data(:, 2);

ids = find(string(bops) == string(bop));
gene = string(genes(ids(1)));
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
legend(h, sprintf('%s (%s): %s', bop, gene, config.gender));
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

