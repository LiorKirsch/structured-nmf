function [W_init, H_init] = get_random_W_H(X, K, parms)
%
% init the random using the seed + the current restart index
   rand_seed = take_from_struct(parms, 'rand_seed', 42);

    rng(rand_seed + parms.restart_ind);
  
  
    W_constraints = take_from_struct(parms, 'W_constraints', 'positive');
    [D, N] = size(X);
    W_init = rand(D, K);
    H_init = rand(K, N) * parms.random_init_spread;
    
    %==== Projection step
    W_init = project_proportions(W_init, W_constraints ,parms);
    Xscale=sum(sum(X));
    Rscale=sum(sum(W_init*H_init));
    sqrnorm=Rscale/Xscale;
    H_init=H_init/sqrnorm;
end

