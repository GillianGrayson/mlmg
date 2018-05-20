clear all;

path = '../../../source/python/D1/gene';

num_points = 20270;

fn = sprintf('%s/spearman_genes_mean.txt', path);
data = importdata(fn);
mean = data.data;

fn = sprintf('%s/spearman_genes_std_from.txt', path);
data = importdata(fn);
std = data.data;

fig = figure;
hLine = scatter(mean(1:num_points), std(1:num_points), 'filled', 'MarkerEdgeColor', 'r', 'SizeData', 5);
title(sprintf('Spearman N=%d', num_points))
xlabel('rho mean', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('rho std', 'Interpreter', 'latex');
box on
propertyeditor(fig)
