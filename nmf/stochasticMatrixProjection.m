function A = stochasticMatrixProjection(A)
% This function projects the matrix A onto the column-stichastic matrices
% It finds the closest l2 matrix that its columns sum to one

    A = projsplx_matrix(A);

%     A_sort = sort(A,1,'descend');
%     for i =1:size(A,2)
%         A(:,i) = projsplx(A(:,i),A_sort(:,i));
%     end


end



function x = projsplx(y,s)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.

m = length(y); bget = false;

% s = sort(y,'descend'); 
tmpsum = 0;

for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - 1)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end;

x = max(y-tmax,0);

end

function A = projsplx_matrix(A)
    % Project a matrix A on the column-stochastic matrices, 
    % such that each column is on the simplex Dn
    % Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
    % 
    % Based on a vector projection implementation by Xiaojing Ye
    % http://www.mathworks.com/matlabcentral/fileexchange/30332-projection-onto-simplex

  
    A_sort = sort(A,1,'descend');

    [num_rows, num_cols] = size(A);
    sorted_cumsum = cumsum(A_sort,1);
    tmax = (sorted_cumsum - ones(size(sorted_cumsum)) ) ./ repmat((1:num_rows)',[1,num_cols]);
    [~,first_val] = max(tmax >= [A_sort(2:end,:);-inf(1,num_cols)]);
    
    select_indcies = num_rows *(0:(num_cols-1)) + first_val;
    zztmax = tmax(select_indcies);


    A = max(A-repmat(zztmax,num_rows,1),zeros(size(A)));

end