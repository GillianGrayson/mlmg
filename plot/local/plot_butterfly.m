function plot_butterfly(plot_data)

% ======== Data for Scatter ========
num_names = size(plot_data.names, 1);
common_names = plot_data.names(plot_data.num_rare + 1:end);
common_metrics_1 = plot_data.metrics_1(plot_data.num_rare + 1:end);
common_metrics_2 = plot_data.metrics_2(plot_data.num_rare + 1:end);
rare_names = plot_data.names(1:plot_data.num_rare);
rare_metrics_1 = plot_data.metrics_1(1:plot_data.num_rare);
rare_metrics_2 = plot_data.metrics_2(1:plot_data.num_rare);

% ======== Data for PDF ========
metrics_1_begin = min(plot_data.metrics_1);
metrics_1_end = max(plot_data.metrics_1);
metrics_1_step = (metrics_1_end - metrics_1_begin) / plot_data.num_bins;
metrics_1_bins = linspace(metrics_1_begin +  0.5 * metrics_1_step, metrics_1_end - 0.5 * metrics_1_step, plot_data.num_bins);

metrics_2_begin = min(plot_data.metrics_2);
metrics_2_end = max(plot_data.metrics_2);
metrics_2_step = (metrics_2_end - metrics_2_begin) / plot_data.num_bins;
metrics_2_bins = linspace(metrics_2_begin +  0.5 * metrics_2_step, metrics_2_end - 0.5 * metrics_2_step, plot_data.num_bins);

metrics_pdf = zeros(plot_data.num_bins);
for id = 1:num_names
    x_id = floor((plot_data.metrics_1(id) - metrics_1_begin) * plot_data.num_bins / (metrics_1_end - metrics_1_begin + 1e-8)) + 1;
    y_id = floor((plot_data.metrics_2(id) - metrics_2_begin) * plot_data.num_bins / (metrics_2_end - metrics_2_begin + 1e-8)) + 1;
    metrics_pdf(x_id, y_id) = metrics_pdf(x_id, y_id) + 1;
end

metrics_pdf = metrics_pdf / (sum(sum(metrics_pdf)) * metrics_1_step * metrics_2_step);
norm = sum(sum(metrics_pdf)) * metrics_1_step * metrics_2_step

% ======== Data for Detla ========
abs_diff = abs(plot_data.metrics_diff);
diff_begin = min(abs_diff);
diff_end = max(abs_diff);
diff_step = (diff_end - diff_begin) / plot_data.num_bins;
diff_bins = linspace(diff_begin +  0.5 * diff_step, diff_end - 0.5 * diff_step, plot_data.num_bins);

diff_pdf = zeros(plot_data.num_bins, 1);
for id = 1:num_names
    x_id = floor((abs_diff(id) - diff_begin) * plot_data.num_bins / (diff_end - diff_begin + 1e-8)) + 1;
    diff_pdf(x_id) = diff_pdf(x_id) + 1;
end
diff_pdf = diff_pdf / (sum(diff_pdf) * diff_step);
norm = sum(diff_pdf) * diff_step

limit_id = floor((abs_diff(plot_data.num_rare) - diff_begin) * plot_data.num_bins / (diff_end - diff_begin + 1e-8)) + 1;

% ======== Scatter figure ========
f1 = figure;
hold all
h = scatter(common_metrics_1, common_metrics_2, 'SizeData', 100, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.5);
set(gca, 'FontSize', 30);
xlabel(plot_data.metrics_1_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(plot_data.metrics_2_label, 'Interpreter', 'latex');

for id = 1:size(rare_names, 1)
    hold all;
    h = scatter(rare_metrics_1(id), rare_metrics_2(id), 'SizeData', 100, 'LineWidth', 5, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.5);
    legend(h, string(rare_names(id)))
end

hold all;
h = plot(xlim, [0 0], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], ylim, 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f1, sprintf('%s/butterfly_scatter_%s.fig', plot_data.save_path, plot_data.suffix))
saveas(f1, sprintf('%s/butterfly_scatter_%s.png', plot_data.save_path, plot_data.suffix))

% ======== PDF figure ========
f2 = figure;
h = imagesc(metrics_1_bins, metrics_2_bins, metrics_pdf');
set(gca, 'FontSize', 30);
xlabel(plot_data.metrics_1_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(plot_data.metrics_2_label, 'Interpreter', 'latex');
colormap hot;
cb = colorbar;
set(gca, 'FontSize', 30);
title(cb, 'PDF');
set(gca,'YDir','normal');

hold all;
h = plot(xlim, [0 0], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], ylim, 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

for id = 1:size(rare_names, 1)
    hold all;
    h = scatter(rare_metrics_1(id), rare_metrics_2(id), 'SizeData', 100, 'LineWidth', 5, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.5);
    legend(h, string(rare_names(id)))
end

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f2, sprintf('%s/butterfly_pdf_%s.fig', plot_data.save_path, plot_data.suffix))
saveas(f2, sprintf('%s/butterfly_pdf_%s.png', plot_data.save_path, plot_data.suffix))

% ======== Delta figure ========
f3 = figure;
hold all
h = plot(diff_bins, diff_pdf, 'LineWidth', 2, 'Color', 'r');
set(gca, 'FontSize', 30);
xlabel('$|\Delta|$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');


hold all;
h = plot(diff_bins(limit_id:end), diff_pdf(limit_id:end), 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(end)], [0 0], 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(limit_id)], [0 diff_pdf(limit_id)], 'LineWidth', 4, 'Color', 'k');

propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f3, sprintf('%s/butterfly_delta_%s.fig', plot_data.save_path, plot_data.suffix))
saveas(f3, sprintf('%s/butterfly_delta_%s.png', plot_data.save_path, plot_data.suffix))

end
