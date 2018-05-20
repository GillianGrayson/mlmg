clear all;

path = '../../../source/python/D1/gene';

num_points = 20270;

fn = sprintf('%s/anova_genes_mean.txt', path);
data = importdata(fn);
mean = data.data;

fn = sprintf('%s/anova_genes_std_from.txt', path);
data = importdata(fn);
std = data.data;

fig = figure;
hLine = scatter(mean(1:num_points), std(1:num_points), 'filled', 'MarkerEdgeColor', 'r', 'SizeData', 5);
title(sprintf('ANOVA N=%d', num_points))
xlabel('pvalue mean', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('pvalue std', 'Interpreter', 'latex');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on
propertyeditor(fig)
