function [new_X, new_W] = add_precision_reg_H(X, W, precision_lambda, ...
    sqrt_precision_matrix,features_mean)

    [num_samp, num_genes] = size(X);
    [num_samp, num_types] = size(W);
    
    mean_matrix = repmat(features_mean, num_types, 1);
    
    new_W_edition = sqrt(precision_lambda) * sqrt_precision_matrix;
    new_X_edition = sqrt_precision_matrix * mean_matrix;

    new_W = [W; new_W_edition];
    new_X = [X ; new_X_edition];
    
end