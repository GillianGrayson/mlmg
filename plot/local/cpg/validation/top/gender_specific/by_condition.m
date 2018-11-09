clear all;

% ======== config ========
config_lvl_1.data_base = 'GSE87571';
config_lvl_1.data_type = 'cpg_data';

config_lvl_1.chromosome_type = 'non_gender';

config_lvl_1.dna_region = 'genic';

config_lvl_1.info_type = 'result';

config_lvl_1.scenario = 'approach';
config_lvl_1.approach = 'top';
config_lvl_1.method = 'linreg_ols';

config_lvl_1.disease = 'any';
config_lvl_1.gender = 'versus';

config_lvl_1.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI')
    config_lvl_1.up = 'D:/YandexDisk/Work/mlmg';
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS')
    config_lvl_1.up = 'D:/Aaron/Bio/mlmg';
else
    config_lvl_1.up = 'E:/YandexDisk/Work/mlmg';
end

% ======== save_config ========
config_lvl_2.data_base = config_lvl_1.data_base;
config_lvl_2.data_type = config_lvl_1.data_type;

config_lvl_2.chromosome_type = config_lvl_1.chromosome_type;

config_lvl_2.dna_region = config_lvl_1.dna_region;

config_lvl_2.info_type = 'result';

config_lvl_2.scenario = 'validation';
config_lvl_2.approach = 'top';
config_lvl_2.method = 'gender_specific';

config_lvl_2.disease = config_lvl_1.disease;
config_lvl_2.gender = 'versus';

config_lvl_2.up = config_lvl_1.up;
config_lvl_2.is_clustering = config_lvl_1.is_clustering;

lvl_1_names = lvl_1_condition(config_lvl_1);
lvl_2_names = lvl_2_condition(config_lvl_1, config_lvl_2);

intersect_names = intersect(lvl_1_names, lvl_2_names);

suffix = sprintf('by_condition');
path = sprintf('%s/data/%s', ...
    config_lvl_1.up, ...
    get_result_path(config_lvl_2));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

d = vertcat('names', intersect_names);

xlswrite(fn, d);






ololo = 1;