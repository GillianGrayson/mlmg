clear all;

% ======== params ========
cpg = 'cg03911306';

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'cpg_data';

config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_ols';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;

config.color = '';

if strcmp(getenv('computername'), 'MSI')
    config.up = 'D:/YandexDisk/Work/mlmg';
else
    config.up = 'E:/YandexDisk/Work/mlmg';
end

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    config.color = 'r';
    plot_linreg_ols_cpg(config, cpg)
    config.gender = 'M';
    config.color = 'b';
    plot_linreg_ols_cpg(config, cpg)
    config.gender = 'versus';
else
    plot_linreg_ols_gene(config, cpg)
end

suffix = sprintf('cpg(%s)', cpg);

up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config));
mkdir(save_path);

box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/linreg_%s.fig', save_path, suffix))
saveas(f, sprintf('%s/linreg_%s.png', save_path, suffix))


function plot_linreg_ols_cpg(config, cpg)

fn = sprintf('%s/data/%s/top.txt', ...
    config.up, ...
    get_result_path(config));

top_data = importdata(fn);

cpgs = string(top_data.textdata);
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
    if mod(num_lines, 1000) == 0
        num_lines = num_lines
    end
    num_lines = num_lines + 1;
end
fclose(fid);

cpg_data_passed = size(indexes, 1);
for id = 1:size(indexes, 1)
    cpg_data_passed(id) = cpg_data(indexes(id));
end

cpg_id = find(cpgs==cpg);

sigma = 3;

slope = slopes(cpg_id);
intercept = intercepts(cpg_id);
intercept_std = intercepts_std(cpg_id);
slope_std = slopes_std(cpg_id);

x_lin = [min(ages), max(ages)];
y_lin = [slope * x_lin(1) + intercept, slope * x_lin(2) + intercept];

slope_minus = slope - sigma * slope_std;
intercept_minus = intercept - sigma * intercept_std;

slope_plus = slope + sigma * slope_std;
intercept_plus = intercept + sigma * intercept_std;

intercept_up = intercept + ((slope_plus * x_lin(2) + intercept_plus) - (slope * x_lin(2) + intercept));
slope_up = slope;

intercept_down = intercept + ((slope_minus * x_lin(2) + intercept_minus) - (slope * x_lin(2) + intercept));
slope_down = slope;

y_up = [slope_up * x_lin(1) + intercept_up, slope_up * x_lin(2) + intercept_up];
y_down = [slope_down * x_lin(1) + intercept_down, slope_down * x_lin(2) + intercept_down];

plot_data.scatter_x = ages_passed;
plot_data.scatter_y = cpg_data_passed;
plot_data.line_x = x_lin;
plot_data.line_y = y_lin;
plot_data.line_y_down = y_down;
plot_data.line_y_up = y_up;
plot_data.line_name = sprintf('%s: %s', cpg, config.gender);
plot_data.color = config.color;

plot_linreg_ols(plot_data)

end

