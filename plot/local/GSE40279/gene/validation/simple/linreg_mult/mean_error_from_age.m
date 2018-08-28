clear all;

base = 'GSE40279';
method = 'enet';
gender_type = 'a';
disease_type = 'any';
gd_approach = 'mean';
gd_validation = 'mean';

fn = sprintf('../../../../data/%s/result/gene/validation/simple/linreg_mult/top/%s/%s/%s/islands_shores/%s/%s/error_from_age.txt', ...
    base, ...
    method, ...
    gender_type, ...
    disease_type, ...
    gd_approach, ...
    gd_validation);
data = importdata(fn);

figure;
h = plot(data(:, 1), data(:, 2));
legend(h, sprintf('gender: %s', gender_type), 'Interpreter', 'latex')
set(gca, 'FontSize', 30);
xlabel('age', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('error', 'Interpreter', 'latex');
hold all;

