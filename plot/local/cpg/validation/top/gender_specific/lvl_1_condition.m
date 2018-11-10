function passed_names = lvl_1_condition(config)
[names, data_1, data_2] = get_gender_specific_data(config);

sigma = 3;

slopes_1 = data_1(:, 3);
slopes_1_std = data_1(:, 5);

slopes_2 = data_2(:, 3);
slopes_2_std = data_2(:, 5);

passed_names = [];
for id = 1:size(names)
    if slopes_1(id) < sigma * slopes_1_std(id) && slopes_2(id) < sigma * slopes_2_std(id)
        passed_names = vertcat(passed_names, names(id));
    end
end

end