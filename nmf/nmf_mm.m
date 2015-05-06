function [W,H,diff_record,time_record]=nmf_mm(parms,X,W_init, H_init)
% The classical NMF which uses multipcative-updates
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% parms : 
%        maxiter : Maximum number of iterations to run
%        loglevel : prints iteration count and changes in connectivity matrix
%        print_interval : print every "print_interval" iterations
%        enforce_prob_norm : If true W is normalized at each iteration
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Based on Lars Kai Hansen, IMM-DTU (c) October 2006
%
% Lior Kirsch 03/2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter = take_from_struct(parms, 'maxiter', 1000);
loglevel = take_from_struct(parms, 'loglevel', 1);
print_interval = take_from_struct(parms, 'print_interval', 100);
W_constraints = take_from_struct(parms, 'W_constraints', 'positive');
record_scores = take_from_struct(parms, 'record_scores', false);
early_stop = take_from_struct(parms, 'early_stop', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(min(X)) < 0
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(X);
W = W_init;
H = H_init;


% use W*H to test for convergence
Xr_old=W*H;

diff_record =nan(1,maxiter);
time_record =nan(1,maxiter); tic;
for iter=1:maxiter

  if (rem(iter,10)==1 && early_stop && (loglevel >0) ) 
    if iter==1,
      gradW = W*(H*H') - X*H';     gradH = (W'*W)*H - W'*X;
      init_grad = norm([gradW; gradH'],'fro');
      init_time=cputime;  
      fprintf('init grad norm %f\n', init_grad);
    end
    if check_for_early_stopping(X,W,H,init_time,init_grad,parms)
        fprintf(' (Iter = %d)\n', iter);
      break
    end
  end

    % Euclidean multiplicative method
    H = H.*(W'*X)./((W'*W)*H+eps);
    W = W.*(X*H')./(W*(H*H')+eps);
    W = project_proportions( W, W_constraints );

    
    % print to screen
    if (rem(iter,print_interval)==0) & (loglevel >0)
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X,W*H);
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
    if record_scores
        diff_record(iter) = nmf_euclidean_dist(X,W*H);
        time_record(iter) = toc;
    end
end
