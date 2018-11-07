clear all;

cpgs = string(importdata('cpgs.txt'));

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

for cpg_id = 1:size(cpgs, 1)
    
    cpg = cpgs(cpg_id)
    
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

end
