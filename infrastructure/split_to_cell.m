function [X_cell, unique_id] = split_to_cell(X, cell_id)
% This function splits the rows of the matrix X in to cells.
%  Input:
%    X
%    cell_id - for each row in X we specify what is the index of the cell.
%
    assert(size(X,1) == length(cell_id), 'number of samples should be the same size(X,1), length(sample_id)' );
    unique_id = unique(cell_id);
    num_cells = length(unique_id);
    X_cell = cell(num_cells,1);
    for i =1:num_cells
       valid_ids = ismember(cell_id, unique_id(i) );
       X_cell{i} = X(valid_ids, :); 
    end
    unique_id = unique_id(:);
end