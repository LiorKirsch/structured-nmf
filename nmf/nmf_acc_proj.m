function [W,H,diff_record,time_record] = nmf_acc_proj(parms,V,Winit,Hinit)
% 
% At each iteration it blocks one matrix solve a least squares solution,
% projects onto the positives and then switches the blocked matrix and does
% the same thing.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% W_init       : N x K (Number of components) 
% H_init       : K x M 
% parms : 
%        maxiter : Maximum number of iterations to run
%        loglevel : prints iteration count and changes in connectivity matrix
%        print_interval : print every "print_interval" iterations
%        enforce_prob_norm : If true W is normalized at each iteration
%        tol: tolerance for a relative stopping condition
%        timelimit, maxiter: limit of time and iterations
%
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Based on "NMF by alternative non-negative least squares using projected
% gradients"    Author: Chih-Jen Lin, National Taiwan University
%
% W,H: output solution

maxiter = take_from_struct(parms, 'maxiter', 1000);
loglevel = take_from_struct(parms, 'loglevel', 1);
print_interval = take_from_struct(parms, 'print_interval', 100);
W_constraints = take_from_struct(parms, 'W_constraints', 'positive');
tol = take_from_struct(parms, 'tolerance', 10^-8);
timelimit = take_from_struct(parms, 'timelimit', 10000); % 2.7 hours
record_scores = take_from_struct(parms, 'record_scores', false);
early_stop = take_from_struct(parms, 'early_stop', true);
accelerated = take_from_struct(parms, 'accelerated', true);

W = Winit; H = Hinit; initt = cputime;

gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;

kkt_dist_H = get_stopping_cond(H,gradH, 'positive');
kkt_dist_W = get_stopping_cond(W',gradW', W_constraints);
initgrad = sqrt( kkt_dist_H^2 + kkt_dist_W ^2 );
% initgrad = norm([gradW; gradH'],'fro');

% fprintf('Init gradient norm %f\n', initgrad); 
tolW = max(0.001,tol)*initgrad; tolH = tolW;

diff_record =nan(1,maxiter);
time_record =nan(1,maxiter); tic;
for iter=1:maxiter,
  % stopping condition
  if early_stop
      
      kkt_dist_H = get_stopping_cond(H,gradH, 'positive');
      kkt_dist_W = get_stopping_cond(W',gradW', W_constraints);
      projnorm = sqrt( kkt_dist_H^2 + kkt_dist_W ^2 );
      
%       projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
      if projnorm < tol*initgrad | cputime-initt > timelimit,
        break;
      end
  end
  
  [H,gradH,iterH] = nlssubprob(V,W,H,tolH,1000,'positive',accelerated);
  if iterH==1,
    tolH = 0.1 * tolH; 
  end

  [W,gradW,iterW] = nlssubprob(V',H',W',tolW,1000,W_constraints,accelerated); W = W'; gradW = gradW';
%   W = project_proportions( W, W_constraints ,parms);

  if iterW==1,
    tolW = 0.1 * tolW;
  end

 

  if (iterW==1 & iterH==1 & tolH + tolW < tol*initgrad),
    fprintf('Failed to move\n'); break;
  end
  if rem(iter,print_interval)==0, fprintf('.'); end
   if record_scores
        diff_record(iter) = nmf_euclidean_dist(V,W*H);
        time_record(iter) = toc;
    end
end

eucl_dist = nmf_euclidean_dist(V,W*H);
fprintf('\nIter = %d Final proj-grad norm %g, ||WH -V||=%g\n', iter, projnorm,eucl_dist);

end

