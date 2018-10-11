function fn = get_gene_data_path(config)
fn = sprintf('%s/%s/%s/%s/%s/%s', ...
    config.data_base, ...
    config.data_type, ...
    config.chromosome_type, ...
    config.geo_type, ...
    config.gene_data_type, ...
    'data');
end