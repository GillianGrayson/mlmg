clear all;

path = '../../../../data/GSE40279/result/gene/validation_linreg';

config = 'enet_plane(mean_der_normed)_islands_shores';

fn = sprintf('%s/%s.txt', path, config);
data = importdata(fn);

num_genes = size(data.data, 1);

vals_main = abs(data.data(:, 1));
vals_aux = abs(data.data(:, 2));

fig = figure;
for g_id = 1:num_genes
    hLine = scatter(vals_main(g_id), vals_aux(g_id), 'Marker', '.', 'SizeData', 500);
    legend(hLine, data.textdata(g_id))
    hold all
end
set(gca, 'FontSize', 30);
xlabel('$|r|$ main', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$|r| aux$', 'Interpreter', 'latex');
box on
propertyeditor(fig)