clear all;

num_genes = 50;

path = 'E:/Work/mlmg/source/python/D1';
fn = sprintf('%s/pvals_genes.txt', path);
data = importdata(fn);
names = data.textdata;
rates = data.data;

[sorted_rates, order] = sort(rates, 'descend');

top_genes = {};
top_rates = [];

for i = 1:num_genes
    if ((rates(order(i+1)) > 1))
        
        top_genes{i} = names{order(i+1)};
        top_rates = vertcat(top_rates, rates(order(i+1)));
        
    end
end
top_genes = top_genes';

fig = figure;
bar(categorical(top_genes), top_rates);

ololo = 0;
