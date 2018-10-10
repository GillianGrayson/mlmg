clear all;

age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

data_base = 'GSE40279';
data_type = 'cpg_data';

chromosome_type = 'non_gender';

dna_region = 'genic';

info_type = 'result';

disease = 'any';
scenario = 'approach';
approach = 'top';
method = 'linreg';

is_clustering = 0;

rare_part = 0.0005;

num_bins = 500;
print_rate = 1000;

metrics_id = 1;
if strcmp(method, 'linreg')
    if is_clustering == 1
        metrics_id = 3;
    else
        metrics_id = 1;
    end
end

up = '../../../../../..';

f_fn = sprintf('%s/data/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/top.txt', ...
    up, ...
    data_base, ...
    data_type, ...
    chromosome_type, ...
    dna_region, ...
    info_type, ...
    disease, ...
    'F', ...
    scenario, ...
    approach, ...
    method);
f_top_data = importdata(f_fn);
f_cpgs = f_top_data.textdata(:, 1);
f_metrics = f_top_data.data(:, metrics_id);

if strcmp(method, 'linreg')
    nothing_to_do = 0;
elseif strcmp(method, 'manova')
    f_metrics = -log10(f_metrics);
end

m_fn = sprintf('%s/data/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/top.txt', ...
    up, ...
    data_base, ...
    data_type, ...
    chromosome_type, ...
    dna_region, ...
    info_type, ...
    disease, ...
    'M', ...
    scenario, ...
    approach, ...
    method);
m_top_data = importdata(m_fn);
m_cpgs = m_top_data.textdata(:, 1);
m_metrics = m_top_data.data(:, metrics_id);

if strcmp(method, 'linreg')
    nothing_to_do = 0;
elseif strcmp(method, 'manova')
    m_metrics = -log10(m_metrics);
end

num_cpgs = size(f_cpgs, 1);

f_metrics_passed = f_metrics;
m_metrics_passed = zeros(num_cpgs, 1);
metrics_diff = zeros(num_cpgs, 1);

for cpg_id = 1:num_cpgs
    f_cpg = f_cpgs(cpg_id);
    m_id = find(m_cpgs==string(f_cpg));
    m_metrics_passed(cpg_id) = m_metrics(m_id);
    metrics_diff(cpg_id) = f_metrics_passed(cpg_id) - m_metrics_passed(cpg_id);
    if mod(cpg_id, print_rate) == 0
        cpg_id = cpg_id
    end
end

[metrics_diff_srt, order] = sort(abs(metrics_diff), 'descend');
cpgs_srt = f_cpgs;
f_metrics_srt = zeros(num_cpgs, 1);
m_metrics_srt = zeros(num_cpgs, 1);
for cpg_id = 1:num_cpgs
    cpgs_srt(cpg_id) = f_cpgs(order(cpg_id));
    f_metrics_srt(cpg_id) = f_metrics_passed(order(cpg_id));
    m_metrics_srt(cpg_id) = m_metrics_passed(order(cpg_id));
end

num_rare = floor(rare_part * num_cpgs);

common_cpgs = cpgs_srt(num_rare + 1:end);
common_f_metrics = f_metrics_srt(num_rare + 1:end);
common_m_metrics = m_metrics_srt(num_rare + 1:end);
common_diff = metrics_diff_srt(num_rare + 1:end);

rare_cpgs = cpgs_srt(1:num_rare);
rare_f_metrics = f_metrics_srt(1:num_rare);
rare_m_metrics = m_metrics_srt(1:num_rare);
rare_diff = metrics_diff_srt(1:num_rare);

figure
hold all;
h = plot([min(f_metrics_srt) max(f_metrics_srt)], [0 0], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(m_metrics_srt) max(m_metrics_srt)], 'LineWidth', 2, 'Color', 'k');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all
h = plot(common_f_metrics, common_m_metrics, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 30);
xlabel('metrics F', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('metrics M', 'Interpreter', 'latex');

for cpg_id = 1:size(rare_cpgs, 1)
    cpg_id = cpg_id
    hold all;
    tareget_cpg = rare_cpgs(cpg_id);
    h = plot(rare_f_metrics(cpg_id), rare_m_metrics(cpg_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, sprintf('%s', string(tareget_cpg)))
end

propertyeditor('on')
box on;

f_metrics_begin = min(f_metrics_srt);
f_metrics_end = max(f_metrics_srt);
f_metrics_step = (f_metrics_end - f_metrics_begin) / num_bins;
f_metrics_bins = linspace(f_metrics_begin +  0.5 * f_metrics_step, f_metrics_end - 0.5 * f_metrics_step, num_bins);

m_metrics_begin = min(m_metrics_srt);
m_metrics_end = max(m_metrics_srt);
m_metrics_step = (m_metrics_end - m_metrics_begin) / num_bins;
m_metrics_bins = linspace(m_metrics_begin +  0.5 * m_metrics_step, m_metrics_end - 0.5 * m_metrics_step, num_bins);

metrics_pdf = zeros(num_bins);
for cpg_id = 1:num_cpgs
    id = floor((f_metrics_srt(cpg_id) - f_metrics_begin) * num_bins / (f_metrics_end - f_metrics_begin + 1e-8)) + 1;
    y_id = floor((m_metrics_srt(cpg_id) - m_metrics_begin) * num_bins / (m_metrics_end - m_metrics_begin + 1e-8)) + 1;
    metrics_pdf(id, y_id) = metrics_pdf(id, y_id) + 1;
end

metrics_pdf = metrics_pdf / (sum(sum(metrics_pdf)) * f_metrics_step * m_metrics_step);
norm = sum(sum(metrics_pdf)) * f_metrics_step * m_metrics_step

figure;
h = imagesc(f_metrics_bins, m_metrics_bins, metrics_pdf');
set(gca, 'FontSize', 30);
xlabel('metrics F', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('metrics M', 'Interpreter', 'latex');
colormap hot;
cb = colorbar;
set(gca, 'FontSize', 30);
title(cb, 'PDF');
set(gca,'YDir','normal');

hold all;
h = plot([min(f_metrics_srt) max(f_metrics_srt)], [0 0], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold all;
h = plot([0 0], [min(m_metrics_srt) max(m_metrics_srt)], 'LineWidth', 2, 'Color', 'g');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

for cpg_id = 1:size(rare_cpgs, 1)
    cpg_id = cpg_id
    hold all;
    tareget_cpg = rare_cpgs(cpg_id);
    h = plot(rare_f_metrics(cpg_id), rare_m_metrics(cpg_id), 'o', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'w');
    legend(h, sprintf('%s', string(tareget_cpg)))
end

propertyeditor('on')

abs_diff = abs(metrics_diff_srt);
diff_begin = min(abs_diff);
diff_end = max(abs_diff);
diff_step = (diff_end - diff_begin) / num_bins;
diff_bins = linspace(diff_begin +  0.5 * diff_step, diff_end - 0.5 * diff_step, num_bins);

diff_pdf = zeros(num_bins, 1);
for cpg_id = 1:num_cpgs
    id = floor((abs_diff(cpg_id) - diff_begin) * num_bins / (diff_end - diff_begin + 1e-8)) + 1;
    diff_pdf(id) = diff_pdf(id) + 1;
end
diff_pdf = diff_pdf / (sum(diff_pdf) * diff_step);
norm = sum(diff_pdf) * diff_step

figure;
hold all
h = plot(diff_bins, diff_pdf, 'LineWidth', 2, 'Color', 'r');
set(gca, 'FontSize', 30);
xlabel('$|\Delta|$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');

limit_id = floor((abs_diff(num_rare) - diff_begin) * num_bins / (diff_end - diff_begin + 1e-8)) + 1;
hold all;
h = plot(diff_bins(limit_id:end), diff_pdf(limit_id:end), 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(end)], [0 0], 'LineWidth', 4, 'Color', 'k');
hold all;
h = plot([diff_bins(limit_id) diff_bins(limit_id)], [0 diff_pdf(limit_id)], 'LineWidth', 4, 'Color', 'k');
propertyeditor('on')


