function [W,H,diff_record,time_record]=nmf_als_with_relations(parms,X_models,relation_matrix_for_H,W_init_model, H_init_models)
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
% Current version supports at most one parent per node (tree structure).
%
assert( parms.W_lambda == 0,'only H structure is supported so far');

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



[N,M]=size(X_all);
[N,K]=size(W_init);

W_models = W_init_model;
H_models = H_init_model;
num_models = size(relation_matrix_for_H,1);
for i=1:num_models
    Xscale=sum(sum(X_models{i}));
    %INIT and rescaling
    Rscale=sum(sum(W_init*H_init));
    sqrnorm=sqrt(Rscale/Xscale);
    H_init=H_init/sqrnorm;
    W_init=W_init/sqrnorm;

    Xr_old = W_init*H_init;

    diff_record =nan(1,maxiter);
    time_record =nan(1,maxiter); tic;


    W_models{i} = W_init{i};
    H_models{i} = H_init{i};
    X_models = cell(num_models,1);
end


for iter=1:maxiter
    for model_iter=1:num_models
    %======
    % for each model solve the problem using the related models as
    % priors
    %======
    
        W = W_models{model_iter};
        H = H_models{model_iter};
        X = X_models{model_iter};
        
    %==== Minization step
        relation_coeff = relation_matrix_for_H(model_iter,:);
        relation_coeff_inds = find(relation_coeff);
        reg_X_for_H = X;
        reg_W_for_H = W;
        for relations_iter = 1:length(relation_coeff_inds)
            relation_ind = relation_coeff_inds(relations_iter);
            [reg_X_for_H, reg_W_for_H] = get_reg_for_H(reg_X_for_H,reg_W_for_H, H_lambda, H_models{relation_ind});
        end
        H = solve_als_for_H(H, reg_W_for_H,reg_X_for_H,als_solver);
        
%         relation_coeff = relation_matrix_for_W(model_iter,:);
%         relation_coeff_inds = find(relation_coeff);
%         reg_X_for_W = X;
%         reg_H_for_W = W;
%         for relations_iter = 1:length(relation_coeff_inds)
%             relation_ind = relation_coeff_inds(relations_iter);
%             [reg_X_for_W, reg_H_for_W] = get_reg_for_H(X',H', W_lambda, W_models{relation_ind});
%             reg_X_for_W = reg_X_for_W';  reg_H_for_W =  reg_H_for_W';
%         end
%         W = solve_als_for_W(W, reg_H_for_W,reg_X_for_W,als_solver);

        W = solve_als_for_W(W, H,X,als_solver);
        %==== Projection step
        W = project_proportions( W, W_constraints,parms );

        
        W_models{model_iter} = W;
        H_models{model_iter} = H;

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
