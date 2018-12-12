clear all;

cpgs = string(importdata('cpgs.txt'));
config.is_plot_regions = 1;

% ======== config ========
config.data_base = 'GSE87571';
config.data_type = 'cpg_data';

config.cross_reactive = 'cross_reactive_excluded';
config.snp = 'snp_excluded';
config.chromosome_type = 'non_gender';

config.dna_region = 'genic';

config.info_type = 'result';

config.scenario = 'approach';
config.approach = 'top';
config.method = 'linreg_variance_ols';

config.disease = 'any';
config.gender = 'versus';

config.is_clustering = 0;

config.color = '';

config.up = get_up_data_path(); 

[annotations_map, labels] = get_annotations(config);

for cpg_id = 1:size(cpgs, 1)
    
    cpg = cpgs(cpg_id)
    
    curr_ann = annotations_map(cpg);
    genes_str = curr_ann(find(labels=="UCSC_REFGENE_NAME"));
    genes = strsplit(genes_str,';');
    genes = unique(genes);
    genes = strjoin(genes, ';');
    
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
    
    suffix = sprintf('cpg(%s)', cpg);
    
    % ======== processing ========
    f = figure;
    if strcmp(config.gender, 'versus')
        config.is_plot_regions = 1;
        config.gender = 'F';
        config.color = 'r';
        plot_linreg_ols_cpg(config, cpg)
        config.gender = 'M';
        config.color = 'b';
        plot_linreg_ols_cpg(config, cpg)
        config.gender = 'versus';
    else
        plot_linreg_ols_cpg(config, cpg)
    end
   
    box on;
    b = gca;
    legend(b,'off');
    yl = ylim;
    if yl(1) < 0
       ylim([0, yl(2)]); 
    end
    if yl(2) > 1
       ylim([yl(1), 1]);
    end
    title(sprintf('%s(%s)', cpg, genes), 'FontSize', 16)
    
    savefig(f, sprintf('%s/%d_linreg_variance_%s.fig', save_path, cpg_id, suffix))
    saveas(f, sprintf('%s/%d_linreg_variance_%s.png', save_path, cpg_id, suffix))
    
    f = figure;
    if strcmp(config.gender, 'versus')
        config.is_plot_regions = 0;
        config.gender = 'F';
        config.color = 'r';
        plot_linreg_variance_ols_cpg(config, cpg)
        config.gender = 'M';
        config.color = 'b';
        plot_linreg_variance_ols_cpg(config, cpg)
        config.gender = 'versus';
    else
        plot_linreg_variance_ols_cpg(config, cpg)
    end
    
    box on;
    b = gca;
    legend(b,'off');
    yl = ylim;
    if yl(1) < 0
       ylim([0, yl(2)]); 
    end
    if yl(2) > 1
       ylim([yl(1), 1]);
    end
    title(sprintf('%s(%s)', cpg, genes), 'FontSize', 16)
    
    savefig(f, sprintf('%s/%d_linreg_variance_diff_%s.fig', save_path, cpg_id, suffix))
    saveas(f, sprintf('%s/%d_linreg_variance_diff_%s.png', save_path, cpg_id, suffix))
    
end
