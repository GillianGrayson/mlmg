clear all;

% ======== params ========
cpg = 'cg22345911';

% ======== config ========
config.data_base = 'GSE40279';
config.data_type = 'cpg_data';

config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;
config.up = '../../../../../..';

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    plot_linreg(config, cpg)
    config.gender = 'M';
    plot_linreg(config, cpg)
    config.gender = 'versus';
else
    plot_linreg(config, cpg)
end

suffix = sprintf('cpg(%s)_gender(%s)', cpg, config.gender);


up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/linreg_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/linreg_%s.png', save_path, suffix))


function plot_linreg(config, cpg)
print_rate = 1000;
fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));

top_data = importdata(fn);
cpgs = top_data.textdata;
slopes = top_data.data(:, 3);
intercepts = top_data.data(:, 4);

indexes = get_attributes_indexes(config);
ages = get_ages(config);

ages_passed = zeros(size(indexes, 1), 1);
for id = 1:size(indexes, 1)
    index = indexes(id);
    ages_passed(id) = ages(index);
end

fn = sprintf('%s/data/%s/average_beta.txt', config.up, config.data_base);
fid = fopen(fn);
num_lines = 1;
while ~feof(fid)
    tline = strsplit(fgetl(fid), '\t');
    curr_cpg = string(tline(1));
    if strcmp(curr_cpg, cpg)
        cpg_data = str2double(tline(2:end))';
        break;
    end
    if mod(num_lines, print_rate) == 0
        num_lines = num_lines
    end
    num_lines = num_lines + 1;
end
fclose(fid);

cpg_data_passed = size(indexes, 1);
for id = 1:size(indexes, 1)
    cpg_data_passed(id) = cpg_data(indexes(id));
end

cpg_id = find(string(cpgs)==string(cpg));

slope = slopes(cpg_id);
intercept = intercepts(cpg_id);
x_lin = [min(ages), max(ages)];
y_lin = [slope * x_lin(1) + intercept, slope * x_lin(2) + intercept];

hold all;
h = plot(ages_passed, cpg_data_passed, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
color = get(h, 'Color');

hold all;
h = plot(x_lin, y_lin, '-', 'LineWidth', 3);
legend(h, sprintf('%s: %s', cpg, config.gender));
set(h, 'Color', color)
set(gca, 'FontSize', 30);
xlabel('age', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\beta$', 'Interpreter', 'latex');

box on;
end

