function cpgs = get_cross_reactive_cpgs(config)

fn = sprintf('%s/data/cross_reactive_cpgs.txt', ...
    config.up);

cpgs = importdata(fn);

end
