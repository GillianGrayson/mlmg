clear all;

% ======== params ========
config.metrics_rank = 1;
config.plot_method = 1;
config.metrics_diff_id = 1;
config.metrics_diff_direction = 'descend';
config.part = 0.05;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'gene_data';

config.chromosome_type = 'non_gender';

config.geo_type = 'islands_shores';
config.gene_data_type = 'mean';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg';

config.disease = 'any';
config.gender = 'any';

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

save_config.geo_type = config.geo_type;
save_config.gene_data_type = config.gene_data_type;

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

% ======== processing =======
gender_specific(config, save_config);