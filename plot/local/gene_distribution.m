clear all;

base = 'GSE40279';
data_type = 'mean';
geo = 'islands_shores';
age_ann = 'age';
genes = {'ELOVL2', ...
    'FHL2', ...
    'TRIM59', ...
    'OTUD7A', ...
    'SOD2', ...
    'TRAPPC3', ...
    'PEX5L', ...
    'USP21'};
genes = genes';

fn = sprintf('../../data/%s/attributes.txt', base);
ann = importdata(fn);
keys = strsplit(string(ann{1}), ' ')';
age_id = 0;
gender_id = 0;
disease_id = 0;
for key_id = 1:size(keys, 1)
    if string(keys{key_id}) == string(age_ann)
        age_id = key_id;
    end
end
ages = zeros(size(ann, 1)-1, 1);
for id = 2:size(ann, 1)
    vals = strsplit(string(ann{id}), ' ')';
    curr_ann = str2double(string(vals{age_id}));
    ages(id-1) = curr_ann;
end

fn = sprintf('../../data/%s/gene_data/%s/%s/gene_data.txt', ...
    base, ...
    data_type, ...
    geo);
data = importdata(fn);
all_genes = data.textdata;
all_data = data.data;

figure;
for gene_id = 1:size(genes, 1)
    gene = string(genes(gene_id));
    idx = find(all_genes==gene);
    gene_data = all_data(idx, :);
    
    h = scatter(ages, gene_data);
    legend(h, gene)
    set(gca, 'FontSize', 30);
    xlabel('age', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('$\beta$', 'Interpreter', 'latex');
    hold all;
end



