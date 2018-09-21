clear all;

age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

base = 'GSE52588';
gender = 'any';
disease = 'any';
geo = 'islands_shores';

fn = sprintf('../../../../../../data/%s/result/cpg/approach/statistics/%s/%s/%s/top.txt', ...
    base, ...
    gender, ...
    disease, ...
    geo);
pdf_data = importdata(fn);

hold all;
h = plot(pdf_data(:, 1), pdf_data(:, 2), 'LineWidth', 3);
legend(h, sprintf('base:%s gender:%s disease:%s', base, gender, disease))
set(gca, 'FontSize', 30);
xlabel('$\beta$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');




