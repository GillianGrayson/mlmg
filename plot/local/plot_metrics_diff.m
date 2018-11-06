function plot_metrics_diff(plot_data)

% ======== Data for Plot ========
num_names = size(plot_data.names, 1);
x = linspace(1, num_names, num_names);
y = plot_data.metrics_diff;

% ======== Plot figure ========
f = figure;

hold all;
h = plot(x, y, '-', 'LineWidth', 3);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(gca, 'FontSize', 30);
xlabel('\#', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('metrics', 'Interpreter', 'latex');
box on;


propertyeditor('on')
box on;
b = gca; legend(b,'off');

savefig(f, sprintf('%s/metrics_diff_%s.fig', plot_data.save_path, plot_data.suffix))
saveas(f, sprintf('%s/metrics_diff_%s.png', plot_data.save_path, plot_data.suffix))

end
