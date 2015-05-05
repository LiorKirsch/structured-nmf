function [W_init, H_init] = get_random_W_H(X,K,parms)

W_constraints = take_from_struct(parms, 'W_constraints', 'positive');

[D,N] = size(X);

W_init = rand(D,K);
H_init = rand(K,N);

%==== Projection step
switch W_constraints
    case 'on_simplex'
        W_init = stochasticMatrixProjection(W_init');
        W_init = W_init';
    case {'inside_simplex', 'positive'}
        % Do nothing since W is randomized between zero and one
    otherwise
        error('unknown option for w_constraints - %s', W_constraints);
end
    

Xscale=sum(sum(X));
Rscale=sum(sum(W_init*H_init));
sqrnorm=Rscale/Xscale;
H_init=H_init/sqrnorm;

end

