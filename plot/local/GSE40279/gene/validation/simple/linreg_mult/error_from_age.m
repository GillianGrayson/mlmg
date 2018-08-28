clear all;

base = 'GSE40279';
method = 'enet';
gender_type = 'any';
disease_type = 'any';
gd_approach = 'mean';
gd_validation = 'mean';

slope = -0.06780602825371794;
intercect = 4.375693445575489;
x_fit = linspace(19,101, 101-19+1);
y_fit = slope * x_fit + intercect;

markers = ['o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h']';
num_markers = size(markers, 1);

fn = sprintf('../../../../data/%s/result/gene/validation/simple/linreg_mult/top/%s/%s/%s/islands_shores/%s/%s/errors.txt', ...
    base, ...
    method, ...
    gender_type, ...
    disease_type, ...
    gd_approach, ...
    gd_validation);
fid = fopen(fn);
tline = fgetl(fid);
num_lines = 1;
figure
h = plot([10, 110], [0 0], 'LineWidth', 0.5, 'Color', 'k');
hold all
h = plot(x_fit, y_fit, 'LineWidth', 2, 'Color', 'r', 'Linestyle', '--');
hold all
title(sprintf('gender: %s', gender_type))
while ischar(tline)
    nums = strsplit(tline)';
    x = str2double(nums(1));
    xs = [];
    ys = [];
    for id = 2:size(nums, 1)
       xs = vertcat(xs, x);
       ys = vertcat(ys, str2double(nums(id))-x);
    end
    
    mk = markers(mod(num_lines, num_markers) + 1);
    
    scatter(xs, ys, 'Marker', mk, 'MarkerFaceColor', 'w');
    set(gca, 'FontSize', 30);
    xlabel('age', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('error', 'Interpreter', 'latex');
    hold all
    
    num_lines = num_lines + 1;
    tline = fgetl(fid);
end



box on


