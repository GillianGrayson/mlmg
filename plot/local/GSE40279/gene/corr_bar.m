clear all;
data = importdata('tmp.txt');
gene_name = 'OTUD7A';

corr = data.data;
name = data.textdata;
x = 1:length(name);

figure
hold on
for i = 1:length(name)
    h = bar(i, corr(i));
    if (mod(i, 2) == 0)
        set(h,'FaceColor','b');
    else
        set(h,'FaceColor','r');
    end
end
hold off


set(gca, 'XTick', x, 'XTickLabel', name);
xtickangle(45)
title(gene_name)
set(gca, 'FontSize', 16);
ylabel('corr coeff', 'Interpreter', 'latex');
hold all;
box on;