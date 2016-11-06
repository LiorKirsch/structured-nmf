function [W,H,diff_record,time_record] = nmf_projections(parms,V,Winit,Hinit)
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
tol = take_from_struct(parms, 'tolerance', 10^-6);
timelimit = take_from_struct(parms, 'timelimit', 10000); % 2.7 hours
record_scores = take_from_struct(parms, 'record_scores', false);
early_stop = take_from_struct(parms, 'early_stop', true);

W = Winit; H = Hinit; initt = cputime;

gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
initgrad = norm([gradW; gradH'],'fro');
fprintf('Init gradient norm %f\n', initgrad); 
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
  
  [H,gradH,iterH] = nlssubprob(V,W,H,tolH,1000,'positive');
  if iterH==1,
    tolH = 0.1 * tolH; 
  end

  [W,gradW,iterW] = nlssubprob(V',H',W',tolW,1000,W_constraints); W = W'; gradW = gradW';
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
fprintf('\nIter = %d Final proj-grad norm %f, ||WH -V||=%f\n', iter, projnorm,eucl_dist);

end

function [H,grad,iter] = nlssubprob(V,W,Hinit,tol,maxiter,W_constraints)
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations

H = Hinit; 
WtV = W'*V;
WtW = W'*W; 

alpha = 1; beta = 0.1;

for iter=1:maxiter,  
  grad = WtW*H - WtV;

  stop_crit = get_stopping_cond(H,grad, W_constraints);
  if stop_crit < tol,
    break
  end

  % projected gradient with step-size search 
  [H, alpha] = proj_grad_with_s_search(H,WtW,grad,alpha, beta,W_constraints);
end

if iter==maxiter,
  fprintf('Max iter in nlssubprob\n');
end
end

function [H, alpha] = proj_grad_with_s_search(H,WtW,grad,alpha, beta,W_constraints)

 Hn = project_proportions( H - alpha*grad, W_constraints);
 
  d = Hn-H;
  gradd=sum(sum(grad.*d));  % grad(f(x_k)) * (x_k+1 - x_k)
  dQd = sum(sum((WtW*d).*d));
  if gradd + 0.5*dQd > 0.01*gradd, 
    % decrease alpha
    while 1,
      alpha = alpha*beta;
      
      Hn = project_proportions( H - alpha*grad, W_constraints);

      d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd <= 0.01*gradd | alpha < 1e-20,      
        H = Hn; break;
      end
    end 
  else 
    % increase alpha
    while 1,
      Hp = Hn;
      alpha = alpha/beta;
      
      Hn = project_proportions( H - alpha*grad, W_constraints);
      
      d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd > 0.01*gradd | Hn == Hp | alpha > 1e10,      
        H = Hp; alpha = alpha*beta; break;
      end
    end 
  end
end

function H = project_proportions( H, projection_type)
% project the proportion onto a set
%   'on_simplex' - projects the W such that each row sums to one
%   'inside_simplex' - projects W such that each row sums is in the interval [0,1]
%   'positive' - make sure each value in W is non-negative
%   'on_simplex_with_noise' - projects the W such that each row sums to one
%        and adds an eps vector such that there would not be a vector which
%        is all zeros

    switch projection_type
        case 'on_simplex'
            H = stochasticMatrixProjection(H,'on');
        case 'on_simplex_with_noise'
            H = stochasticMatrixProjection(H,'on');
            H = H + rand(size(H))* eps; % moves it away from zero
        case 'inside_simplex'
            H = stochasticMatrixProjection(H,'inside');
        case 'positive'
            H=(H>0).*H;
        otherwise
            error('unknown option for projection type - %s', projection_type);
    end

end

