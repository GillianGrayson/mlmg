clear all;

base = 'GSE40279';
method = 'linreg';
gender_type = 'F';
disease_type = 'any';
gd_approach = 'mean';
geo = 'islands_shores';

cluster_id = 1;
num_genes = 500;
x = linspace(1, num_genes, num_genes);

markers = ['o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h']';
num_markers = size(markers, 1);
colors = ['y', 'm', 'c', 'r', 'g', 'b', 'k']';
num_colors = size(colors, 1);

fn = sprintf('../../../../../../data/%s/result/gene/approach/top/%s/%s/%s/%s/%s/top.txt', ...
    base, ...
    method, ...
    gender_type, ...
    disease_type, ...
    gd_approach, ...
    geo);

data = importdata(fn);

metrics_data = data.data(1:num_genes, 3);

cluster_data = data.data(1:num_genes, cluster_id);
max_cluster = max(cluster_data);
min_cluster = min(cluster_data);
num_clusters = max_cluster - min_cluster + 1;

curr_cluster = cluster_data(1);
curr_sorted = 0;
sorted_clusters = zeros(size(cluster_data, 1), 1);
for i = 1:size(cluster_data, 1)
    if(cluster_data(i) ~= curr_cluster)
        curr_cluster = cluster_data(i);
        curr_sorted = curr_sorted + 1;
    end
    sorted_clusters(i) = curr_sorted;
end

figure
subplot(2, 1, 1)
h = plot(x, metrics_data, 'Marker', 'o', 'MarkerFaceColor', 'w');
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
title(sprintf('%s', method))
set(gca, 'FontSize', 30);
xlabel('n', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('metrics', 'Interpreter', 'latex');
hold all

subplot(2, 1, 2)
for c_id = 1:num_clusters
    cluster_curr = min_cluster + (c_id - 1);
    xs = [];
    ys = [];
    for p_id = 1:num_genes
        x_curr = x(p_id);
        y_curr = sorted_clusters(p_id);
        if y_curr == cluster_curr
            xs = vertcat(xs, x_curr);
            ys = vertcat(ys, y_curr);
        end
    end
    
    marker_curr = markers(mod(cluster_curr, num_markers) + 1);
    color_curr = colors(mod(cluster_curr, num_colors) + 1);
    
    h = scatter(xs, ys, 'Marker', marker_curr, 'MarkerEdgeColor', color_curr, 'MarkerFaceColor', 'w');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    set(gca, 'FontSize', 30);
    xlabel('n', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('cluster', 'Interpreter', 'latex');
    hold all
end
box on








