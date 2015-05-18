function corr_matrix = callibrated_corr(matrix_a, matrix_b, corr_type)

assert(size(matrix_a,1) == size(matrix_b,1),'matrix A and B should have same size first dim');
corr_matrix = nan(size(matrix_a,2), size(matrix_b,2));

for i = 1:size(matrix_a,2)
    for j = 1:size(matrix_b,2)
        vec_a = matrix_a(:,i);
        vec_b = matrix_b(:,j);
        vec_b = match_expression_across_experiments(vec_a, vec_b);
        corr_matrix(i,j) = corr(vec_a, vec_b, 'type',corr_type);
    end
end



end