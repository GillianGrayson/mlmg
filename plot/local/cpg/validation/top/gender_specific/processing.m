clear all;

% ======== params ========
config.metrics_rank = 1;
config.plot_method = 1;
config.metrics_diff_id = 2;
config.metrics_diff_direction = 'ascend';
config.part = 0.0005;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'cpg_data';

config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_ols_wo_outliers';

config.disease = 'any';
config.gender = 'versus';

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

save_config.dna_region = config.dna_region;

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
gender_specific(config, save_config);