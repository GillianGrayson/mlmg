function diff_metrics = get_gender_specific_diff_metrics(config)

diff_metrics = zeros(size(config.names, 1), 1);
for gene_id = 1:size(config.names, 1)
    diff_metrics(gene_id) = config.f_metrics(gene_id) - config.m_metrics(gene_id);
end

end