clear all;

base = 'GSE40279';
method = 'enet';
gender_type = 'F';
disease_type = 'any';
gd_approach = 'mean';
gd_validation = 'mean';

tops = 5:5:350;
tops = tops';

corrs = zeros(size(tops, 1), 1);
for top_id = 1:size(tops, 1)
    
    top = tops(top_id);
    fn = sprintf('../../../../data/%s/result/gene/validation/simple/linreg_mult/top/%s/%s/%s/islands_shores/%s/%s/metrics_%d.txt', ...
        base, ...
        method, ...
        gender_type, ...
        disease_type, ...
        gd_approach, ...
        gd_validation, ...
        top);
    data = importdata(fn);
    corrs(top_id) = data.data(1);
end


figure;
h = plot(tops, corrs);
title(sprintf('gender: %s', gender_type), 'Interpreter', 'latex');
legend(h, base)
set(gca, 'FontSize', 30);
xlabel('number of genes', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('corr', 'Interpreter', 'latex');
hold all;

