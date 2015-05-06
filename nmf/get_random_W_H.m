function [W_init, H_init] = get_random_W_H(X,K,parms)

W_constraints = take_from_struct(parms, 'W_constraints', 'positive');

[D,N] = size(X);

W_init = rand(D,K);
H_init = rand(K,N);

%==== Projection step
W_init = project_proportions( W_init, W_constraints ,parms);
    

Xscale=sum(sum(X));
Rscale=sum(sum(W_init*H_init));
sqrnorm=Rscale/Xscale;
H_init=H_init/sqrnorm;

end

