function [W,H,diff_record,time_record]=nmf_als(parms, X, W_init, H_init)
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
%
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Lior Kirsch 03/2015

maxiter = take_from_struct(parms, 'maxiter', 1000);
miniter = take_from_struct(parms, 'miniter', 10);
tol = take_from_struct(parms, 'tolerance', 10^-6);
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

H_markers = take_from_struct(parms, 'H_markers', false(size(H_init)) );

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

assert( all(size(H_init) == size(H_markers)) , 'H_markers should be a boolean array with same size as H');
assert( all(sum(H_markers,1) <=1)  ,'A marker should only be present for a single type');
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

gradW = W*(H*H') - X*H'; gradH = (W'*W)*H - W'*X;
initgrad = norm([gradW; gradH'],'fro');
tolW = max(0.001,tol)*initgrad; tolH = tolW;

all_marker_indices = any(H_markers,1);

for iter=1:maxiter

      % ==== Early stopping is based on distance from the closest x which 
    % ==== 	satisfies the KKT conditions
    if (rem(iter,50)==1 && early_stop) 
        
          [reg_X_for_H, reg_W_for_H] = get_reg_for_H(X(:,~all_marker_indices),W, H_lambda, H_prior(:,~all_marker_indices));
          gradH = (reg_W_for_H'*reg_W_for_H)*H - W'*reg_X_for_H;
          kkt_dist_H = get_stopping_cond(H,gradH, 'positive');
          
          [reg_X_for_W, reg_H_for_W] = get_reg_for_H(X',H', W_lambda, W_prior);
          reg_X_for_W = reg_X_for_W';  reg_H_for_W =  reg_H_for_W';
          gradW = W*(reg_H_for_W*reg_H_for_W') - reg_X_for_W*reg_H_for_W'; 
          kkt_dist_W = get_stopping_cond(W',gradW', W_constraints);
          
          kkt_dist = sqrt( kkt_dist_H^2 + kkt_dist_W ^2 );
        if iter==1,  
          init_kkt_dist = kkt_dist;
        end
    
    
        if (kkt_dist < tol*init_kkt_dist)
            fprintf('Early stopping:  (Iter = %d) (kkt-dist = %g) (HL=%4.2g)', iter, kkt_dist,parms.H_lambda);
            fprintf(' (||X-WH||_F %g)\n', nmf_euclidean_dist(X,W*H) );
            break;
        end
    end
%==== Minization step

    % split H to markers and non markers
 
    H_marker_part = update_H_markers(X(:,all_marker_indices), W, H_markers(:,all_marker_indices),H_lambda,H_prior(:,all_marker_indices));
    [reg_X_for_H, reg_W_for_H] = get_reg_for_H(X(:,~all_marker_indices),W, H_lambda, H_prior(:,~all_marker_indices));
    [H_non_marker_parts, iterH] = solve_als_for_H(H, reg_W_for_H,reg_X_for_H,als_solver,tolH,maxiter,'positive');
    % join H-markers part and H-non-markers part
    H = nan(size(H_init));
    H(:,all_marker_indices) = H_marker_part;
    H(:,~all_marker_indices) = H_non_marker_parts;
    
    [reg_X_for_W, reg_H_for_W] = get_reg_for_H(X',H', W_lambda, W_prior);
    reg_X_for_W = reg_X_for_W';  reg_H_for_W =  reg_H_for_W';
    [W,iterW] = solve_als_for_W(W, reg_H_for_W,reg_X_for_W,als_solver,tolW,maxiter,W_constraints);

%==== Projection step
    W = project_proportions( W, W_constraints ,parms);
   
    
    if iterH <= miniter
        tolH = 0.1 * tolH; 
    end
    if iterW <= miniter
        tolW = 0.1 * tolW;
    end



    % print to screen
    if (rem(iter,print_interval)==0) && (loglevel >0)
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        eucl_dist  = nmf_euclidean_dist(X,W*H);
        eucl_dist_old  = nmf_euclidean_dist(X,Xr_old);
        errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
        fprintf('Iter %d ||X-WH||=%g, mean(X-WH)=%g, diff(k,k-1)=%g, ||X-WH||delta=%g\n',...
            iter,eucl_dist,errorx, diff, eucl_dist_old - eucl_dist);
         
        Xr_old = Xr;

%         eucl_dist  = nmf_euclidean_dist(X,W*H);
%         errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
%         fprintf('Iter %d ||X-WH||=%g, mean(X-WH)=%g, diff(k,k-1)=%g\n',...
%             iter,eucl_dist,errorx, diff);
        if errorx < 10^(-5), break, end
    end

     if record_scores
        diff_record(iter) = nmf_euclidean_dist(X,W*H);
        time_record(iter) = toc;
     end
    
  
end

if (iter==maxiter)
    disp(sprintf('max limit iteration (%d) kkt-dist %g reached (HL=%4.2g)', ...
                 maxiter, kkt_dist, parms.H_lambda));
end
end

%%%%=== add parms
function [H,num_iters] = solve_als_for_H(init_H, W,X,als_solver, tol,maxiter,W_constraints)
    switch als_solver
        case 'pinv_project'
            H=(W*pinv(W'*W))'*X;
            H=H.*(H>0);
            num_iters = nan;
        case 'active_set'
            [H,gradHX,num_iters] = nnlsm_activeset(W,X,1,0,init_H);
        case 'blockpivot'
            [H,gradHX,num_iters] = nnlsm_blockpivot(W,X,0,init_H);
        case 'accel_proj'
            [H,grad,num_iters] = nlssubprob_accel(X,W,init_H,tol,maxiter,W_constraints);
        case 'cjlin'
            accelerated = false;
            [H,grad,num_iters] = nlssubprob(X,W,init_H,tol,maxiter,W_constraints,accelerated);
        otherwise 
            error('unknown solver - %s', als_solver);
    end
end
function [W,num_iters] = solve_als_for_W(init_W, H,X,als_solver, tol,maxiter,W_constraints)
    switch als_solver
        case 'pinv_project'
            W = ((pinv(H*H')*H)*X')';
            num_iters = nan;
%===== I removed this because this is done in a next when I project to or on the simplex.
%             W=(W>0).*W;  
%=====
        case 'active_set'
            [W,gradW,num_iters] = nnlsm_activeset(H',X',1,0,init_W'); 
            W=W'; 
        case 'blockpivot'
            [W,gradW,num_iters] = nnlsm_blockpivot(H',X',0,init_W');
            W=W';
        case 'accel_proj'
            [W,grad,num_iters] = nlssubprob_accel(X',H',init_W',tol,maxiter,W_constraints);
            W=W';
        case 'cjlin'
            accelerated = false;
            [W,grad,num_iters] = nlssubprob(X',H',init_W',tol,maxiter,W_constraints,accelerated);
            W=W';
        otherwise 
            error('unknown solver - %s', als_solver);
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
