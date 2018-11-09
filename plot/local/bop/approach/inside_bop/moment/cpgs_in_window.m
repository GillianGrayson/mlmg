clear all;

% ======== params ========
bop = 'chr1:118727816-118728097*S_Shore';
start_cpg = 'cg09497789';

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

config.color = '';

if strcmp(getenv('computername'), 'MSI') 
    config.up = 'D:/YandexDisk/Work/mlmg'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    config.up = 'D:/Aaron/Bio/mlmg'; 
else 
    config.up = 'E:/YandexDisk/Work/mlmg'; 
end

% ======== processing ========
f = figure;
if strcmp(config.gender, 'versus')
    config.gender = 'F';
    config.color = 'r';
    plot_cpg_in_window(config, bop, start_cpg)
    config.gender = 'M';
    config.color = 'b';
    plot_cpg_in_window(config, bop, start_cpg)
    config.gender = 'versus';
else
    plot_linreg_gene(config, cpg)
end

bop_name = strrep(bop, ':', '_');
bop_name = strrep(bop_name, '*', '_');
bop_name = strrep(bop_name, '-', '_');
suffix = sprintf('bop(%s)_start_cpg(%s)', bop_name, start_cpg);

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

box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/%s.fig', save_path, suffix))
saveas(f, sprintf('%s/%s.png', save_path, suffix))


function plot_cpg_in_window(config, bop, start_cpg)

fn = sprintf('%s/data/%s/inside_bop.txt', ...
    config.up, ...
    get_result_path(config));

inside_bop_data = importdata(fn, ' ');
bops = inside_bop_data.textdata(:, 1);
cpgs = inside_bop_data.textdata(:, 3);

ids = find(string(bops) == string(bop));
num_cpgs_in_bop = size(ids, 1);
bop_cpgs = {};
target_cpgs = {};
for i = 1:num_cpgs_in_bop
    bop_cpgs{i} = string(cpgs(ids(i)));
    if strcmp(bop_cpgs{i}, start_cpg)
        target_cpgs{1} = string(cpgs(ids(i)));
        target_cpgs{2} = string(cpgs(ids(i+1)));
        target_cpgs{3} = string(cpgs(ids(i+2)));
    end
end

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
cpg_data = {};
cpg_order = zeros;
cpgs_passed = 0;
while ~feof(fid)
    tline = strsplit(fgetl(fid), '\t');
    curr_cpg = string(tline(1));
    if strcmp(curr_cpg, target_cpgs{1}) || strcmp(curr_cpg, target_cpgs{2}) || strcmp(curr_cpg, target_cpgs{3})
        cpg_order(cpgs_passed+1) = find([target_cpgs{:}] == curr_cpg);
        cpg_data{cpgs_passed+1} = str2double(tline(2:end))';
        cpgs_passed = cpgs_passed + 1;
        if cpgs_passed == 3
            break;
        end
    end
    if mod(num_lines, 1000) == 0
        num_lines = num_lines
    end
    num_lines = num_lines + 1;
end
fclose(fid);

for cpg_id = 1:3

    cpg_data_passed = size(indexes, 1);
    for id = 1:size(indexes, 1)
        cpg_data_passed(id) = cpg_data{cpg_order(cpg_id)}(indexes(id));
    end

    plot_data.scatter_x = ages_passed;
    plot_data.scatter_y = cpg_data_passed;
    plot_data.color = config.color;

    hold all;
    subplot(3, 1, cpg_id);
    h = scatter(plot_data.scatter_x, plot_data.scatter_y, 'MarkerEdgeColor', plot_data.color, 'SizeData', 100, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.5);
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    set(gca, 'FontSize', 24);
    xlabel('age', 'Interpreter', 'latex');
    set(gca, 'FontSize', 24);
    ylabel('$\beta$', 'Interpreter', 'latex');
    xlim([min(ages)-1, max(ages)+1])
    title(target_cpgs{cpg_order(cpg_id)});

    box on;
end

end

