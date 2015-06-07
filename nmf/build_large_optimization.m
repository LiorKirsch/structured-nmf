function [X_out, W_out] = build_large_optimization(X_cell, W_cell, relation_matrix, h_lambda)

    % each X_cell    is   samp X genes
    % each W_cell    is   samp X types
    % relation_matrix   is  num_regions X num_regions
    
    num_regions = length(X_cell);
    [num_samp,num_genes] = size(X_cell{1});
    [num_samp,num_types] = size(W_cell{1});
    
    X_out = [];
    W_out = [];
    
    for i = 1:num_regions
       X_out = cat(2, X_out, X_cell{i}); 
       W_out = cat(2, W_out, W_cell{i}); 
    end

    for i = 1:num_regions
        for j = i+1:num_regions
            multp = h_lambda * relation_matrix(i,j);
            new_W = zeros(num_types, num_types * num_regions);
            new_W( :, 1 + (i-1)*num_types) = multp *eye(num_types);    
            new_W( :, 1 + (j-1)*num_types) = -multp *eye(num_types);    
            
            new_X = zeros(num_types, num_genes * num_regions);

            X_out = cat(1, X_out, new_X); 
            W_out = cat(1, W_out, new_W); 
        end
    end



end