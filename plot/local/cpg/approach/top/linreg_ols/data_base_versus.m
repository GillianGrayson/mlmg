clear all;

cpgs = string(importdata('cpgs.txt'));
data_bases = ["GSE40279"; "GSE87571"];

for cpg_id = 1:size(cpgs, 1)
    
    cpg = cpgs(cpg_id)
    
    f = figure;
    for data_base_id = 1:size(data_bases, 1)
        
        % ======== config ========
        config.data_base = data_bases(data_base_id);
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
        subplot(1, size(data_bases, 1), data_base_id)
        
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
        
        box on;
        b = gca;
        legend(b,'off');
        title(config.data_base, 'FontSize', 16)
        yl = ylim;
        if yl(1) < 0
            ylim([0, yl(2)]);
        end
        if yl(2) > 1
            ylim([yl(1), 1]);
        end
        sgtitle(cpg, 'FontSize', 20)
        propertyeditor('on')
        
        suffix = sprintf('cpg(%s)_data_bases(%s)', cpg, join(sort(data_bases), '_'));
        
        config.data_base = "versus";
        up_save = 'C:/Users/user/Google Drive/mlmg/figures';
        
        save_path = sprintf('%s/%s', ...
            up_save, ...
            get_result_path(config));
        mkdir(save_path);
        
        savefig(f, sprintf('%s/linreg_ols_%s.fig', save_path, suffix))
        saveas(f, sprintf('%s/linreg_ols_%s.png', save_path, suffix))
        
    end
    
end