function [W,H,diff_record,time_record]=nmf_als(parms,X,W_init, H_init)
% At each iteration it blocks one matrix solve a least squares solution,
% projects onto the positives and then switches the blocked matrix and does
% the same thing.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% W_init       : N x K (Number of components) 
% H_init       : K x M 
% W_lambda     : (from parms)
% H_lambda     : (from parms)
% W_prior      : N x k (from parms  1 <= k <= K)
% H_prior      : k x M (from parms  1 <= k <= K)
% parms : 
%        maxiter : Maximum number of iterations to run
%        loglevel : prints iteration count and changes in connectivity matrix
%        print_interval : print every "print_interval" iterations
%        enforce_prob_norm : If true W is normalized at each iteration
%
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Based on Lars Kai Hansen, IMM-DTU (c) October 2006
%
% Lior Kirsch 03/2015

maxiter = take_from_struct(parms, 'maxiter', 1000);
loglevel = take_from_struct(parms, 'loglevel', 1);
print_interval = take_from_struct(parms, 'print_interval', 100);
W_constraints = take_from_struct(parms, 'W_constraints', 'positive');
record_scores = take_from_struct(parms, 'record_scores', false);
early_stop = take_from_struct(parms, 'early_stop', true);
als_solver = take_from_struct(parms, 'als_solver', 'blockpivot');

W_lambda = take_from_struct(parms, 'W_lambda', 0);
H_lambda = take_from_struct(parms, 'H_lambda', 0);
W_prior = take_from_struct(parms, 'W_prior', nan);
H_prior = take_from_struct(parms, 'H_prior', nan);

if (W_lambda>0)
    if ( isnan(W_prior) )
        W_prior = zeros(size(W_init));
        disp('W_prior not specified using zero regularization as default');
    end
    assert( size(W_prior,2) == size(W_init,2) , 'W_prior should have the same number of instances as W (dim=1)');
else
    W_prior = nan(size(W_init));
end
if (H_lambda>0)
    if ( isnan(H_prior) )
        H_prior = zeros(size(H_init));
        disp('H_prior not specified using zero regularization as default');
    end
    assert( size(H_prior,2) == size(H_init,2) , 'H_prior should have the same number of features as H (dim=2)');
else
    H_prior = nan(size(H_init));   
end

assert( ~(W_lambda>0 && H_lambda>0), 'priors for both H and W is not supported');

[N,M]=size(X);
[N,K]=size(W_init);
Xscale=sum(sum(X));
%INIT
W = W_init;
H = H_init;
Rscale=sum(sum(W*H));
sqrnorm=sqrt(Rscale/Xscale);
H=H/sqrnorm;
W=W/sqrnorm;

Xr_old = W*H;

diff_record =nan(1,maxiter);
time_record =nan(1,maxiter); tic;



for iter=1:maxiter

    if (rem(iter,10)==1 && early_stop) 
        if iter==1,
          gradW = W*(H*H') - X*H';     gradH = (W'*W)*H - W'*X;
          init_grad = norm([gradW; gradH'],'fro');
          init_time=cputime;  
%           fprintf('init grad norm %f\n', init_grad);
        end
        if check_for_early_stopping(X,W,H,init_time,init_grad,parms)
            fprintf(' (Iter = %d)\n', iter);
          break
        end
    end

%==== Minization step
    [reg_X_for_H, reg_W_for_H] = get_reg_for_H(X,W, H_lambda, H_prior);
    H = solve_als_for_H(H, reg_W_for_H,reg_X_for_H,als_solver);
    [reg_X_for_W, reg_H_for_W] = get_reg_for_H(X',H', W_lambda, W_prior);
    reg_X_for_W = reg_X_for_W';  reg_H_for_W =  reg_H_for_W';
    W = solve_als_for_W(W, reg_H_for_W,reg_X_for_W,als_solver);

%==== Projection step
    switch W_constraints
        case 'on_simplex'
            W = stochasticMatrixProjection(W');
            W = W';
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
            error('unknown option for w_constraints - %s', W_constraints);
    end

    % print to screen
    if (rem(iter,print_interval)==0) && (loglevel >0)
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist  = nmf_euclidean_dist(X,W*H);
        errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
        fprintf('Iter %d ||X-WH||=%g, mean(X-WH)=%g, diff(k,k-1)=%g\n',...
            iter,eucl_dist,errorx, diff);
        if errorx < 10^(-5), break, end
    end

     if record_scores
        diff_record(iter) = nmf_euclidean_dist(X,W*H);
        time_record(iter) = toc;
    end
end

if (iter==maxiter)
    disp('max limit iteration reached');
end
end


function H = solve_als_for_H(init_H, W,X,als_solver)
    switch als_solver
        case 'pinv_project'
            H=(W*pinv(W'*W))'*X;
            H=H.*(H>0);
        case 'active_set'
            [H,gradHX,subIterH] = nnlsm_activeset(W,X,1,0,init_H);
        case 'blockpivot'
            [H,gradHX,subIterH] = nnlsm_blockpivot(W,X,0,init_H);
    end
end
function W = solve_als_for_W(init_W, H,X,als_solver)
    switch als_solver
        case 'pinv_project'
            W = ((pinv(H*H')*H)*X')';
%===== I removed this because this is done in a next when I project to or on the simplex.
%             W=(W>0).*W;  
%=====
        case 'active_set'
            [W,gradW,subIterW] = nnlsm_activeset(H',X',1,0,init_W'); 
            W=W'; 
        case 'blockpivot'
            [W,gradW,subIterW] = nnlsm_blockpivot(H',X',0,init_W');
            W=W'; 
    end
end

function [reg_X_for_H, reg_W_for_H] = get_reg_for_H(X,W,lambda_H, H_regularizer)
    k = size(W,2);
    k_prior = size(H_regularizer,1);
    if (lambda_H > 0)
        reg_W_for_H = [W;lambda_H*eye(k_prior,k)];
        reg_X_for_H = [X;lambda_H*H_regularizer];
    else
        reg_W_for_H = W;
        reg_X_for_H = X;
    end
               
end
