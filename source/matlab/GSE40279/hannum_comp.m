clear all;

num_top = 66;

path = '../../python/Methylation/cpg';
fn = sprintf('%s/enet_bootstrap_cpgs.txt', path);
data = importdata(fn);
cpgs = string(data.textdata(1:end, 1));
cpgs = cpgs(1:num_top);
genes_raw = string(data.textdata(1:end, 2));
genes_raw = genes_raw(1:num_top);
genes = [];
for i = 1:size(genes_raw, 1)
    tmp = string(strsplit(genes_raw(i), ';'));
    for j = 1:size(tmp, 2)
        genes = vertcat(genes, tmp(j));
    end
end
genes = unique(genes, 'stable');

path = '../../../data/GSE40279';
fn = sprintf('%s/hannum_best_cpgs.txt', path);
data = importdata(fn);
cpgs_hannum = string(data);
cpgs_int = intersect(cpgs, cpgs_hannum, 'stable')
num_int_hannum = size(cpgs_int, 1)

path = '../../../data/GSE40279';
fn = sprintf('%s/table.txt', path);
data = importdata(fn);
genes_claudio = string(data);
genes_int = intersect(genes, genes_claudio, 'stable');
num_int_claudio = size(genes_int, 1)
