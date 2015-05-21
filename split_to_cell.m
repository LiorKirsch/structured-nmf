function [X_cell, unique_id] = split_to_cell(X, sample_id)
% 
    assert(size(X,1) == length(sample_id), 'number of samples should be the same size(X,1), length(sample_id)' );
    unique_id = unique(sample_id);
    num_cells = length(unique_id);
    X_cell = cell(num_cells,1);
    for i =1:num_cells
       valid_ids = ismember(sample_id, unique_id(i) );
       X_cell{i} = X(valid_ids, :); 
    end
    unique_id = unique_id(:);
end