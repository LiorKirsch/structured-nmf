function W = project_proportions( W, projection_type )
% project the proportion onto a set
%   'on_simplex' - projects the W such that each row sums to one
%   'inside_simplex' - projects W such that each row sums is in the interval [0,1]
%   'positive' - make sure each value in W is non-negative
%   'on_simplex_with_noise' - projects the W such that each row sums to one
%        and adds an eps vector such that there would not be a vector which
%        is all zeros

    switch projection_type
        case 'on_simplex'
            W = stochasticMatrixProjection(W');
            W = W';
        case 'on_simplex_with_noise'
            W = stochasticMatrixProjection(W');
            W = W';
            W = W + rand(size(W))* eps;
        case 'inside_simplex'
            W=(W>0).*W;
            row_sum = sum(W,2);
            valid_rows =  (0 <= row_sum) & (row_sum <= 1);
            W_tmp = stochasticMatrixProjection(W(~valid_rows,:)');
            W(~valid_rows,:) = W_tmp';
        case 'positive'
             if strcmp(als_solver, 'pinv_project')
                W=(W>0).*W;
            end
        otherwise
            error('unknown option for w_constrains - %s', projection_type);
    end

end

