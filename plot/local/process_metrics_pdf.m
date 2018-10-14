function pdf = process_metrics_pdf(target_pdf, config)
pdf = target_pdf;
if strcmp(config.method, 'manova')
    pdf = log10(target_pdf + 1.0e-16);
end
end
