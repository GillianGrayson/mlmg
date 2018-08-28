clear all;

age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

base = 'GSE40279';
gender = 'F';
disease_type = 'healthy';

age_lim = 55;

num_cpgss = 5;

fn = sprintf('../../../data/%s/result/cpg/approach/bend/linreg/%s/any/bend_%d.txt', ...
    base, ...
    gender, ...
    age_lim);
bend_data = importdata(fn);

cpgs = bend_data.textdata(1:num_cpgss, 1);

fn = sprintf('../../../data/%s/attributes.txt', base);
ann = importdata(fn);

keys = strsplit(string(ann{1}), ' ')';
age_id = 0;
gender_id = 0;
disease_id = 0;
for key_id = 1:size(keys, 1)
    if string(keys{key_id}) == string(age_ann)
        age_id = key_id;
    end
    if string(keys{key_id}) == string(gender_ann)
        gender_id = key_id;
    end
    if string(keys{key_id}) == string(disease_ann)
        disease_id = key_id;
    end
end

indexes = [1:size(ann, 1)-1]';

if gender_id > 0  
    genders = [];
    for id = 2:size(ann, 1)
        vals = strsplit(string(ann{id}), ' ')';
        curr_ann = string(vals{gender_id});
        genders = vertcat(genders, curr_ann);
    end
    
    if ~strcmp(string(gender), 'any')
        curr_indexes = [];
        for id = 1:size(genders, 1)
            if strcmp(genders(id), gender)
                curr_indexes = vertcat(curr_indexes, id);
            end
        end
        tmp_indexes = intersect(indexes, curr_indexes);
        indexes = tmp_indexes;
    end
end

if disease_id > 0  
    diseases = [];
    for id = 2:size(ann, 1)
        vals = strsplit(string(ann{id}), ' ')';
        curr_ann = string(vals{disease_id});
        diseases = vertcat(diseases, curr_ann);
    end
    
    if ~strcmp(string(disease_type), 'any')
        curr_indexes = [];
        for id = 1:size(diseases, 1)
            if strcmp(diseases(id), disease_type)
                curr_indexes = vertcat(curr_indexes, id);
            end
        end
        tmp_indexes = intersect(indexes, curr_indexes);
        indexes = tmp_indexes;
    end
end

ages = zeros(size(ann, 1)-1, 1);
for id = 2:size(ann, 1)
    vals = strsplit(string(ann{id}), ' ')';
    curr_ann = str2double(string(vals{age_id}));
    ages(id-1) = curr_ann;
end

indexes_less = [];
for id = 1:size(indexes, 1)
   index = indexes(id);
   if (ages(index) <  age_lim)
       indexes_less = vertcat(indexes_less, index);
   end
end

indexes_more = [];
for id = 1:size(indexes, 1)
   index = indexes(id);
   if (ages(index) >=  age_lim)
       indexes_more = vertcat(indexes_more, index);
   end
end

ages_less = zeros(size(indexes_less, 1), 1);
for id = 1:size(indexes_less, 1)
   index = indexes_less(id);
   ages_less(id) = ages(index);
end

ages_more = zeros(size(indexes_more, 1), 1);
for id = 1:size(indexes_more, 1)
   index = indexes_more(id);
   ages_more(id) = ages(index);
end

fn = sprintf('../../../data/%s/result/cpg/approach/bend/linreg/%s/any/bend_data_%d.txt', ...
    base, ...
    gender, ...
    age_lim);
data = importdata(fn);
cpgs_names = data.textdata;
cpgs_data = data.data;

figure;
colors = {};
for cpg_id = 1:size(cpgs, 1)
    cpg_name = string(cpgs(cpg_id));
    idx = find(cpgs_names==cpg_name);
    cpg_data = cpgs_data(idx, :)';
    
    cpg_data_less = size(indexes_less, 1);
    for id = 1:size(indexes_less, 1)
        cpg_data_less(id) = cpg_data(indexes_less(id));
    end
    
    cpg_data_more = size(indexes_more, 1);
    for id = 1:size(indexes_more, 1)
        cpg_data_more(id) = cpg_data(indexes_more(id));
    end
    
    slope_less = bend_data.data(cpg_id, 2);
    intercept_less = bend_data.data(cpg_id, 3);
    x_lin_less = [min(ages_less), max(ages_less)];
    y_lin_less = [slope_less * x_lin_less(1) + intercept_less, slope_less * x_lin_less(2) + intercept_less];
    
    slope_more = bend_data.data(cpg_id, 7);
    intercept_more = bend_data.data(cpg_id, 8);
    x_lin_more = [min(ages_more), max(ages_more)];
    y_lin_more = [slope_more * x_lin_more(1) + intercept_more, slope_more * x_lin_more(2) + intercept_more];
    
    hold all;
    h = plot(ages_less, cpg_data_less, 'o');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    
    color = get(h, 'Color');
    colors = vertcat(colors, color);
    
    hold all;
    h = plot(ages_more, cpg_data_more, 'x', 'Color', color);
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    
    hold all;
    h = plot(x_lin_less, y_lin_less, '-', 'LineWidth', 3, 'Color', color);
    legend(h, cpg_name);
    
    hold all;
    h = plot(x_lin_more, y_lin_more, '-', 'LineWidth', 3, 'Color', color);
     set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    
    set(gca, 'FontSize', 30);
    xlabel('age', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('$\beta$', 'Interpreter', 'latex');

end

box on;




