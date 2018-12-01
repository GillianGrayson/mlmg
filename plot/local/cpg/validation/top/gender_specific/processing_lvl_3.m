clear all;

num_lvl_1 = 1;
num_lvl_2 = 1;
target_method = "linreg_ols";

data_base = 'GSE87571';
data_type = 'cpg_data';

cross_reactive = 'cross_reactive_excluded';
snp = 'snp_excluded';
chromosome_type = 'non_gender';

dna_region = 'genic';

info_type = 'result';

disease = 'any';
gender = 'versus';

lvl_1_scenario = 'approach';
lvl_1_approach = 'top';
lvl_1_methods = ["linreg_ols"];
lvl_1_suffixes = ["_outliers_limit(0.3)_outliers_sigma_(2.0)"];
lvl_1_experiments = [5];

lvl_2_scenario = 'validation';
lvl_2_approach = 'top';
lvl_2_methods = ["gender_specific"];
lvl_2_suffixes = [""];
lvl_2_experiments = [5];

all_metrics_labels = [];
intersection_names = [];
lvl_1_metrics_map = {};
lvl_2_metrics_map = {};

for lvl_1_id = 1:num_lvl_1
    
    clearvars config;
    
    config.data_base = data_base;
    config.data_type = data_type;
    
    config.cross_reactive = cross_reactive;
    config.snp = snp;
    config.chromosome_type = chromosome_type;
    
    config.dna_region = dna_region;
    
    config.info_type = info_type;
    
    config.scenario = lvl_1_scenario;
    config.approach = lvl_1_approach;
    config.method = lvl_1_methods(lvl_1_id);
    
    config.disease = disease;
    config.gender = gender;
    
    config.is_clustering = 0;
    
    config.up = get_up_data_path();
    
    config.suffix = lvl_1_suffixes(lvl_1_id);
    
    config.experiment = lvl_1_experiments(lvl_1_id);
    
    [names, metrics_labels, metrics_map] = lvl_1_condition(config);
    
    if lvl_1_id == 1 
        intersection_names = names;
    else
        intersection_names = intersect(intersection_names, names);
    end
    
    all_metrics_labels = horzcat(all_metrics_labels, metrics_labels);
    lvl_1_metrics_map{end + 1} = metrics_map;
        
end

for lvl_2_id = 1:num_lvl_2
    
    clearvars config;
    
    config.data_base = data_base;
    config.data_type = data_type;
    
    config.cross_reactive = cross_reactive;
    config.snp = snp;
    config.chromosome_type = chromosome_type;
    
    config.dna_region = dna_region;
    
    config.info_type = info_type;
    
    config.scenario = lvl_2_scenario;
    config.approach = lvl_2_approach;
    config.method = lvl_2_methods(lvl_2_id);
    
    config.disease = disease;
    config.gender = gender;
    
    config.is_clustering = 0;
    
    config.up = get_up_data_path();
    
    config.suffix = lvl_2_suffixes(lvl_2_id);
    
    config.experiment = lvl_2_experiments(lvl_2_id);
    
    [names, metrics_labels, metrics_map] = lvl_2_condition(target_method, config);
    
    intersection_names = intersect(intersection_names, names);
    
    all_metrics_labels = horzcat(all_metrics_labels, metrics_labels);
    lvl_2_metrics_map{end + 1} = metrics_map;
    
    save_config = config;
        
end

metrics_data = [];
for name_id = 1:size(intersection_names,1)
    name = string(intersection_names(name_id));
    data = [];
    for lvl_1_id = 1:num_lvl_1
        data = horzcat(data, lvl_1_metrics_map{lvl_1_id}(name));
    end
    
    for lvl_2_id = 1:num_lvl_2
        data = horzcat(data, lvl_2_metrics_map{lvl_2_id}(name));
    end
    metrics_data = vertcat(metrics_data, data);
end

lvl_1_suffix = sprintf('%d', lvl_1_experiments(1));
for lvl_1_id = 2:num_lvl_1
    lvl_1_suffix = sprintf('%s_%d', lvl_1_suffix, lvl_1_experiments(lvl_1_id));
end

lvl_2_suffix = sprintf('%d', lvl_2_experiments(1));
for lvl_2_id = 2:num_lvl_2
    lvl_2_suffix = sprintf('%s_%d', lvl_2_suffix, lvl_2_experiments(lvl_2_id));
end

suffix = sprintf('lvl_3_by_lvl_1(%s)_lvl_2(%s)', lvl_1_suffix, lvl_2_suffix);
path = sprintf('%s/data/%s', ...
    save_config.up, ...
    get_result_path(save_config));
mkdir(path)
fn = sprintf('%s/%s.xlsx', ...
    path, ...
    suffix);

d = vertcat("names", intersection_names);
for metrics_id = 1:size(metrics_data, 2)
    d = horzcat(d, vertcat(all_metrics_labels(metrics_id), string(metrics_data(:, metrics_id))));
end

xlswrite(fn, d);