function [real_scores, rand_scores] = results_struct_to_mat(results)
%
% flatten the recursive results structure into a 3-way matrix
% scores has dims:  [num_lambdas, num_regions, num_types]
% 
    [num_lambdas, num_regions, num_types] = get_dims(results);
    real_scores = zeros(num_lambdas, num_regions, num_types);
    rand_scores = zeros(num_lambdas, num_regions, num_types);
    i_var_outer = 1; 
    results_outer = results{i_var_outer}; 
        
    for i_lambda = 1: num_lambdas
        r = results_outer{i_lambda};        
        for i_region = 1:num_regions
            real_scores(i_lambda, i_region, 1:num_types) = r.celltype_scores{i_region};
            rand_scores(i_lambda, i_region, 1:num_types) = ...
                r.randbase_celltype_score{i_region};
        end
    end
end

% ==========================================
function [num_lambdas, num_regions, num_types] = get_dims(results)
    num_lambdas = length(results{1});
    r = results{1}{1};
    num_regions = length(r.regions);
    num_types = length(r.celltype_scores{1}); 
end



