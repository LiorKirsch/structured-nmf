function [W_init, H_init] = get_random_W_H(X, K, parms)
%
% init the random using the seed + the current restart index
   rand_seed = take_from_struct(parms, 'rand_seed', 42);
   init_by = take_from_struct(parms, 'init_by', 'random');

    rng(rand_seed + parms.restart_ind);
  
  
    W_constraints = take_from_struct(parms, 'W_constraints', 'positive');
    [D, N] = size(X);
    W_init = rand(D, K);
   
    
    %==== Projection step
    W_init = project_proportions(W_init, W_constraints ,parms);
    switch init_by
        case 'random'
            H_init = rand(K, N);
            Xscale=sum(sum(X));
            Rscale=sum(sum(W_init*H_init));
            sqrnorm=Rscale/Xscale;
            H_init=H_init/sqrnorm;
        case 'samples'
            selected_inds = randperm(D,K);
            H_init = X(selected_inds, :);
        case 'samples_with_noise'
            selected_inds = randperm(D,K);
            H_init = X(selected_inds, :);
            H_init = H_init .* rand(K, N);
        case 'combination'
            combine_matrix = rand(K, D);
            combine_matrix = project_proportions(combine_matrix, 'on_simplex_with_noise' ,parms);
            H_init = combine_matrix * X;
        otherwise
            error('unkown init method - %s', init_by);
    end
    H_init = H_init * parms.random_init_spread;
end

