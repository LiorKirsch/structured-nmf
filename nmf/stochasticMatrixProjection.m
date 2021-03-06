function A = stochasticMatrixProjection(A,projection_type)
% This function projects the matrix A onto the column-stochastic matrices.
% It finds the closest matrix that its columns sum to one, and each 
% entry in the matrix is in the interval [0,1]
%
% input:
%     A - matrix, which columns to project
%     projection_type - which type of projection to make:
%       'on'     - on the simplex, each column sums to one (default)
%       'inside' - inside the simplex, each column sum is in the interval[0,1]
%
% example:
%     vec_2d = rand(2,200)*3 - 1;
%     proj_vec = stochasticMatrixProjection(vec_2d,'on');
%     subplot(1,2,1); 
%     plot([vec_2d(1,:);proj_vec(1,:)],[vec_2d(2,:);proj_vec(2,:)],... 
%          [vec_2d(1,:);proj_vec(1,:)],[vec_2d(2,:);proj_vec(2,:)],'.'); axis equal
%     subplot(1,2,2); 
%     proj_vec = stochasticMatrixProjection(vec_2d,'inside');
%     plot([vec_2d(1,:);proj_vec(1,:)],[vec_2d(2,:);proj_vec(2,:)],... 
%          [vec_2d(1,:);proj_vec(1,:)],[vec_2d(2,:);proj_vec(2,:)],'.'); axis equal
%
% 
%  Lior Kirsch   May, 2015.
%

if ~exist('projection_type','var')
    projection_type = 'on';
end

switch projection_type
    case 'on'
        A = projsplx_matrix(A);
    case 'inside'
        % First project the entries to non-negative.
        % Then, for rows which are still not inside the simplex do the
        % 'on-simplex' projection.
        
        A=(A>0).*A;
        col_sum = sum(A,1);
        valid_cols =  (0 <= col_sum) & (col_sum <= 1);
        A_tmp = projsplx_matrix(A(:,~valid_cols));
        A(:,~valid_cols) = A_tmp;
    otherwise
        error('projection type is either "on" or "inside"');
end

end


function A = projsplx_matrix(A)
    % Project a matrix A on the column-stochastic matrices, 
    % such that each column is on the simplex Dn
    % Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
    % 
    % Based on a vector projection implementation by Xiaojing Ye
    % http://www.mathworks.com/matlabcentral/fileexchange/30332-projection-onto-simplex
    % Algorithm is explained as in http://arxiv.org/abs/1101.6081

  
    A_sort = sort(A,1,'descend');

    [num_rows, num_cols] = size(A);
    sorted_cumsum = cumsum(A_sort,1);
%     tmax = (sorted_cumsum - ones(size(sorted_cumsum)) ) ./ repmat((1:num_rows)',[1,num_cols]);
    tmax = bsxfun(@rdivide, sorted_cumsum -1, (1:num_rows)');
    
    [~,first_val] = max(tmax >= [A_sort(2:end,:);-inf(1,num_cols)]);
    
    select_indcies = num_rows *(0:(num_cols-1)) + first_val;
    zztmax = tmax(select_indcies);

    
%     A = max(A-repmat(zztmax,num_rows,1),zeros(size(A)));
    A = max(bsxfun(@minus, A, zztmax),0); 

end