clear all;

base = 'GSE87571';
age_ann = 'age';
gender_ann = 'gender';
disease_ann = 'disease';

gender_type = 'M';
disease_type = 'any';

edges = 0:5:110;

fn = sprintf('../../../../data/%s/attributes.txt', base);
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
    
    if ~strcmp(string(gender_type), 'any')
        curr_indexes = [];
        for id = 1:size(genders, 1)
            if strcmp(genders(id), gender_type)
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

ages_passed = zeros(size(indexes, 1), 1);
for id = 1:size(indexes, 1)
   index = indexes(id);
   ages_passed(id) = ages(index);
end

figure;
h = histogram(ages_passed, edges);
title(sprintf('gender: %s', gender_type), 'Interpreter', 'latex');
legend(h, base)
set(gca, 'FontSize', 30);
xlabel('$age$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$count$', 'Interpreter', 'latex');
hold all;

