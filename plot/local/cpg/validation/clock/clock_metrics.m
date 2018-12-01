clear all;

metric = "r_test";

data_base = "GSE87571";
data_type = "cpg_data";

cross_reactive = "cross_reactive_excluded";
snp = "snp_excluded";
chromosome_type = "non_gender";

dna_region = "genic";

info_type = "result";

disease = "any";
gender = "versus";

lvl_1_scenario = "approach";
lvl_1_approach = "top";
lvl_1_method = "custom";
lvl_1_suffix = "";

lvl_2_scenario = "validation";
lvl_2_approach = "clock";
lvl_2_method = "linreg_mult";
lvl_2_suffix = "";

% ======== config_lvl_1 ========
config_lvl_1.data_base = data_base;
config_lvl_1.data_type = data_type;

config_lvl_1.cross_reactive = cross_reactive;
config_lvl_1.snp = snp;
config_lvl_1.chromosome_type = chromosome_type;

config_lvl_1.dna_region = dna_region;

config_lvl_1.info_type = info_type;

config_lvl_1.scenario = lvl_1_scenario;
config_lvl_1.approach = lvl_1_approach;
config_lvl_1.method = lvl_1_method;

config_lvl_1.disease = disease;
config_lvl_1.gender = gender;

config_lvl_1.is_clustering = 0;

config_lvl_1.up = get_up_data_path();

config_lvl_1.suffix = lvl_1_suffix;

% ======== config_lvl_2 ========
config_lvl_2.data_base = data_base;
config_lvl_2.data_type = data_type;

config_lvl_2.cross_reactive = cross_reactive;
config_lvl_2.snp = snp;
config_lvl_2.chromosome_type = chromosome_type;

config_lvl_2.dna_region = dna_region;

config_lvl_2.info_type = info_type;

config_lvl_2.scenario = lvl_2_scenario;
config_lvl_2.approach = lvl_2_approach;
config_lvl_2.method = lvl_2_method;

config_lvl_2.disease = disease;
config_lvl_2.gender = gender;

config_lvl_2.is_clustering = 0;

config_lvl_2.up = get_up_data_path();

config_lvl_2.suffix = lvl_2_suffix;


% ======== processing ========
f = figure;
if strcmp(config_lvl_2.gender, 'versus')
    config_lvl_2.gender = 'F';
    config_lvl_2.color = 'r';
    plot_clock_metrics(config_lvl_1, config_lvl_2, metric)
    config_lvl_2.gender = 'M';
    config_lvl_2.color = 'b';
    plot_clock_metrics(config_lvl_1, config_lvl_2, metric)
    config_lvl_2.gender = 'versus';
else
    plot_clock_metrics(config_lvl_1, config_lvl_2, metric)
end

suffix = sprintf('metric(%s)', metric);

up_save = 'C:/Users/user/Google Drive/mlmg/figures';

save_path = sprintf('%s/%s', ...
    up_save, ...
    get_result_path(config_lvl_2));
mkdir(save_path);

box on;
b = gca;
legend(b,'off');
xlim([1 100])
title(sprintf('clock_method(%s)', config_lvl_1.method), 'FontSize', 16, 'Interpreter', 'none')

savefig(f, sprintf('%s/clock_method(%s)_metric(%s).fig', save_path, config_lvl_1.method, metric))
saveas(f, sprintf('%s/clock_method(%s)_metric(%s).png', save_path, config_lvl_1.method, metric))