clear all;

path = '../../../source/python/D1/gene';

num_points = 1000;

fn = sprintf('%s/linreg_genes_mean.txt', path);
data = importdata(fn);
mean = data.data(:, 1);

fn = sprintf('%s/linreg_genes_std_from.txt', path);
data = importdata(fn);
std = data.data(:, 1);

fig = figure;
hLine = scatter(mean(1:num_points), std(1:num_points), 'filled', 'MarkerEdgeColor', 'r', 'SizeData', 5);
title(sprintf('LinReg N=%d', num_points))
xlabel('rho mean', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('rho std', 'Interpreter', 'latex');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on
propertyeditor(fig)
