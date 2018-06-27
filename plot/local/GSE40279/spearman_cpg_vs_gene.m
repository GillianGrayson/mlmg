clear all;

path = '../../../source/python/D1/cpg_gene';

num_points = 15000;

fn = sprintf('%s/spearman_spec.txt', path);
data = importdata(fn);
mean = data.data;

fn = sprintf('%s/spearman_spec_from_mean.txt', path);
data = importdata(fn);
mean_gene = data.data;

fn = sprintf('%s/spearman_spec_from_std.txt', path);
data = importdata(fn);
std_gene = data.data;

fig = figure;
hLine = scatter(mean(1:num_points), mean_gene(1:num_points), 'filled', 'MarkerEdgeColor', 'r', 'SizeData', 5);
title(sprintf('Spearman N=%d', num_points))
xlabel('rho cpg', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('rho mean gene', 'Interpreter', 'latex');
box on
propertyeditor(fig)

fig = figure;
hLine = scatter(mean(1:num_points), std_gene(1:num_points), 'filled', 'MarkerEdgeColor', 'r', 'SizeData', 5);
title(sprintf('Spearman N=%d', num_points))
xlabel('rho cpg', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('rho std gene', 'Interpreter', 'latex');
box on
propertyeditor(fig)