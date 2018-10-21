function diff_metrics = get_gender_specific_diff_metrics(config)

diff_metrics = zeros(size(config.names, 1), 1);
for gene_id = 1:size(config.names, 1)
    diff_metrics(gene_id) = config.metrics_1(gene_id) - config.metrics_2(gene_id);
end

end