function diff = process_diff(a, b, diff_type)
diff = a - b;
if diff_type == 1
    diff = a - b;
elseif diff_type == 2
    diff = sqrt(a^2 + b^2);
end
end
