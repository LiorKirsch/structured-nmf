function [W_init, H_init] = get_random_W_H(X,K,parms)

enforce_W_prob_norm = take_from_struct(parms, 'enforce_W_prob_norm', false);

[D,N] = size(X);

W_init = rand(D,K);
H_init = rand(K,N);


Xscale=sum(sum(X));
if enforce_W_prob_norm
    W_init=W_init./(repmat(sum(W_init,2),1,size(W_init,2))+eps); % normalize rows to unit length
    Rscale=sum(sum(H_init));
    sqrnorm=sqrt(D*Rscale/Xscale);
    H_init=H_init/sqrnorm;
else
    Rscale=sum(sum(W_init*H_init));
    sqrnorm=sqrt(Rscale/Xscale);
    H_init=H_init/sqrnorm;
    W_init=W_init/sqrnorm;
end





end

