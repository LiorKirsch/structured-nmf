function [new_X, new_W] = add_self_reg_H(X, W, old_H, self_lambda, ...
    self_matrix,self_reg_type, features_mean, features_std)
% Won't really work...
% I need some way of normalizing it.
% 
% even if I had a way of normalizng it it would still be strange since
% K (num_types) is usually very small , so it is mean, std, corr over a
% very small number of samples
% 
% *  I can divide each gene by the var (not std) of H_old
% *  What happens when the std( H_old(gene_i) )  =  0 ???
%

    [num_types, num_genes] = size(old_H);
    
    
    switch self_reg_type
        case 'pre_calc_mean'
            old_H_centered = bsxfun(@minus, old_H, features_mean);
            old_H_normalized = bsxfun(@rdivide, old_H_centered, features_std.^2);
            
            new_W_edition = sqrt(self_lambda) * old_H_normalized;
            new_X_edition = self_matrix + bsxfun(@times, old_H_normalized, features_mean);

            new_W = [W; new_W_edition];
            new_X = [X ; new_X_edition];
        case 'data_mean'
            mean_transform = eye(num_types) - 1/num_types * ones(num_types);
            old_H_centered = bsxfun(@minus, old_H, mean(old_H));
            new_edition = sqrt(self_lambda) * mean_transform *old_H_centered;
            new_edition = bsxfun(@rdivide, new_edition, features_std.^2);
            new_W = [W; new_edition];
            new_X = [X ; self_matrix];
        otherwise
            error('unkown type of self regulation - %s', self_reg_type);
            
    end
            

end