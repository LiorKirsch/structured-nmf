function show_proportions(proportions, cell_types)
    [num_types,num_samples] = size(proportions);
    assert(length(cell_types) == num_types,'should be the same number of cell types');
    for i =1:num_types
        figure('Name',sprintf('%s proportions',cell_types{i}) );
        hist(proportions(i,:), 10);
        title(cell_types{i} )
    end
end