clear all;

path = '../../../../data/GSE40279/result/gene/validation_linreg';

config = 'enet_metrics(mean)_vals(mean)';

fn = sprintf('%s/%s.txt', path, config);
data = importdata(fn);

num_genes = size(data.data, 1);

corr_coeff = data.data(:, 2);
genes_count = linspace(1, num_genes, num_genes);

fig = figure;
hLine = plot(genes_count, abs(corr_coeff))
xlabel('#', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$|r|$', 'Interpreter', 'latex');
box on
propertyeditor(fig)