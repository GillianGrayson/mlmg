clear all;

path = '../../../data/GSE40279';

suffix = 'mean_shores_n';

fn = sprintf('%s/R2s_%s.txt', path, suffix);
data = importdata(fn);

nums = data(:, 1);
r2s = data(:, 2);

fig = figure;
hLine = plot(nums, r2s);
legend(hLine, suffix)
xlabel('num', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$R^2$', 'Interpreter', 'latex');
box on
propertyeditor(fig)