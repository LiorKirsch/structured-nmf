function W = project_proportions( W, projection_type ,parms)
% project the proportion onto a set
%   'on_simplex' - projects the W such that each row sums to one
%   'inside_simplex' - projects W such that each row sums is in the interval [0,1]
%   'positive' - make sure each value in W is non-negative
%   'on_simplex_with_noise' - projects the W such that each row sums to one
%        and adds an eps vector such that there would not be a vector which
%        is all zeros

    switch projection_type
        case 'on_simplex'
            W = stochasticMatrixProjection(W','on');
            W = W';
        case 'on_simplex_with_noise'
            W = stochasticMatrixProjection(W','on');
            W = W';
            W = W + rand(size(W))* eps; % moves it away from zero
        case 'inside_simplex'
            W = stochasticMatrixProjection(W','inside');
            W = W';
        case 'positive'
            % other methods (other then 'pinv_project' take care of the
            % non-negativity so there is no need to project again here
            if strcmp(parms.nmf_method, 'pinv_project')
                W=(W>0).*W;
            end
        otherwise
            error('unknown option for projection type - %s', projection_type);
    end

end

