clear all;

path = '../../../../../../data/GSE40279/result/gene/validation/simple/match/mean/islands';
fn = 'match_matrix.txt';
fn = sprintf('%s/%s', path, fn);
data = importdata(fn);

a = data{1};
