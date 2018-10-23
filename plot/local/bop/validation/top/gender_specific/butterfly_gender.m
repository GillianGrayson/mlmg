clear all;

% ======== params ========
config.metrics_rank = 1;
config.plot_method = 2;
config.part = 0.005;

config.num_bins = 200;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'bop_data';

config.chromosome_type = 'non_gender';

config.class_type = 'ClassAB';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'manova';

config.disease = 'any';
config.gender = '';

config.is_clustering = 0;

if strcmp(getenv('computername'), 'MSI')
    config.up = 'D:/YandexDisk/Work/mlmg';
else
    config.up = 'E:/YandexDisk/Work/mlmg';
end

% ======== save_config ========
save_config.data_base = config.data_base;
save_config.data_type = config.data_type;

save_config.chromosome_type = config.chromosome_type;

save_config.class_type = config.class_type;

save_config.info_type = 'result';

save_config.scenario = 'validation';
save_config.approach = 'top';
save_config.method = 'gender_specific';

save_config.disease = config.disease;
save_config.gender = 'versus';

save_config.up = 'C:/Users/user/Google Drive/mlmg/figures';
save_config.is_clustering = config.is_clustering;

% ======== processing ========
butterfly(config, save_config);