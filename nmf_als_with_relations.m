function [W_models,H_models,diff_record_models,time_record]=nmf_als_with_relations(parms,X_models,relation_matrix_for_H,W_init_model, H_init_model)
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

H_lambda_priors = take_from_struct(parms, 'H_lambda_prior', zeros(length(H_init_model),1));
H_prior_models = take_from_struct(parms, 'H_prior_models', repmat({nan},length(H_init_model),1) );

H_markers_models = take_from_struct(parms, 'H_markers', ...
    cellfun(@(x) false(size(x)), H_init_model ,'UniformOutput',false) );


% if (W_lambda>0)
%     if ( isnan(W_prior) )
%         W_prior = zeros(size(W_init));
%         disp('W_prior not specified using zero regularization as default');
%     end
%     assert( size(W_prior,2) == size(W_init,2) , 'W_prior should have the same number of instances as W (dim=1)');
% else
%     W_prior = cellfun(@(x) nan(size(x)),W_init,'UniformOutput',false);
% end
if any(H_lambda_priors>0)
    if ~iscell(H_prior_models)
        disp('Using the same priors for all regions');
        H_prior_models = repmat({H_prior_models}, length(H_init_model),1);
    end
    if length(H_lambda_priors) == 1
       disp('Using the same priors for all regions');
       H_lambda_priors = H_lambda_priors *  ones(length(H_init_model),1);
    end
    for i =1:length(H_init_model)
        assert( all(size(H_init_model{i}) == size(H_prior_models{i})) ,...
            'H_model_prior (%d) should have the same dimentions as H',i);
    end

end

assert( all(size(H_init_model{1}) == size(H_markers_models{1})) , 'H_markers should be a boolean array with same size as H');
assert( all(sum(H_markers_models{1},1) <=1)  ,'A marker should only be present for a single type');
assert( ~(W_lambda>0 && H_lambda>0), 'priors for both H and W is not supported');



[N,M]=cellfun(@size,X_models);
% [N,K]=size(W_init);

W_models = W_init_model;
H_models = H_init_model;
num_models = size(relation_matrix_for_H,1);

assert(num_models == length(X_models),'the number of elements in X and in the relation matrix are not the same');
for i=1:num_models
    W_init = W_models{i};
    H_init = H_models{i};
    
    Xscale=sum(sum(X_models{i}));
    
    %INIT and rescaling
    Rscale=sum(sum(W_init*H_init));
    sqrnorm=sqrt(Rscale/Xscale);
    H_init=H_init/sqrnorm;
    W_init=W_init/sqrnorm;

%     Xr_old = W_init*H_init;

    time_record =nan(1,maxiter); tic;

    W_models{i} = W_init;
    H_models{i} = H_init;
end

Xr_old_models = cellfun(@(W,H) W*H,W_models,H_models,'UniformOutput',false);
diff_record_models = repmat( {nan(1,maxiter)}, num_models,1);
init_grad = cell(num_models,1);
init_time = cell(num_models,1);

for iter=1:maxiter
    
     if (rem(iter,10)==1 && early_stop) 
        reached_early = false(num_models,1);
        for model_iter=1:num_models
            W = W_models{model_iter};
            H = H_models{model_iter};
            X = X_models{model_iter};
            
            if iter==1,
              gradW = W*(H*H') - X*H';     gradH = (W'*W)*H - W'*X;
              init_grad{model_iter} = norm([gradW; gradH'],'fro');
              init_time{model_iter} = cputime;  
    %           fprintf('init grad norm %f\n', init_grad);
            end
            reached_early(model_iter) =  check_for_early_stopping(X,W,H,init_time{model_iter},init_grad{model_iter},parms);
        end
        
        if all(reached_early)
            fprintf(' (Iter = %d)\n', iter);
          break
        end
     end
    
     
    for model_iter=1:num_models
    %======
    % for each model solve the problem using the related models as
    % priors
    %======
    
        W = W_models{model_iter};
        H = H_models{model_iter};
        X = X_models{model_iter};
        H_prior = H_prior_models{model_iter};
        H_markers = H_markers_models{model_iter};
        H_lambda_prior = H_lambda_priors(model_iter);
        Xr_old = Xr_old_models{model_iter};
        diff_record = diff_record_models{model_iter};
    %==== Minization step
    
        relation_coeff = relation_matrix_for_H(model_iter,:);
        relation_coeff_inds = find(relation_coeff);
       
        
        % split H to markers and non markers
        all_marker_indices = any(H_markers,1);
        
        % First, update the marker genes
        tmp_priors = zeros(size(H(:,all_marker_indices)));
        sum_regulerizer = sum(relation_coeff(relation_coeff_inds));
        for relations_iter = 1:length(relation_coeff_inds)
            relation_ind = relation_coeff_inds(relations_iter);
            curr_regulerizer = relation_coeff(relation_ind) / sum_regulerizer;
            tmp_priors = tmp_priors + curr_regulerizer*H_models{relation_ind}(:,all_marker_indices) ;
        end
        H_marker_part = update_H_markers(X(:,all_marker_indices), W, H_markers(:,all_marker_indices),H_lambda*sum_regulerizer,tmp_priors);
        
        reg_X_for_H = X(:,~all_marker_indices);
        reg_W_for_H = W;
        % Then, update the rest of the genes
        for relations_iter = 1:length(relation_coeff_inds)
            relation_ind = relation_coeff_inds(relations_iter);
            curr_regulerizer = H_lambda * relation_coeff(relation_ind);
            [reg_X_for_H, reg_W_for_H] = get_reg_for_H(reg_X_for_H,...
                reg_W_for_H, curr_regulerizer, H_models{relation_ind}(:,~all_marker_indices));
        end
        
        % add the priors
        if H_lambda_prior > 0;
            [reg_X_for_H, reg_W_for_H] = get_reg_for_H(reg_X_for_H,...
                    reg_W_for_H, H_lambda_prior, H_prior(:,~all_marker_indices));
        end 
        H_non_marker_parts = solve_als_for_H(H, reg_W_for_H,reg_X_for_H,als_solver);
    
        % join H-markers part and H-non-markers part
        H = nan(size(H_init));
        H(:,all_marker_indices) = H_marker_part;
        H(:,~all_marker_indices) = H_non_marker_parts;
    

        W = solve_als_for_W(W, H,X,als_solver);
        %==== Projection step
        W = project_proportions( W, W_constraints,parms );

        
        W_models{model_iter} = W;
        H_models{model_iter} = H;

        % print to screen
        if (rem(iter,print_interval)==0) && (loglevel >0)
            Xr = W*H;
            diff = sum(sum(abs(Xr_old-Xr)));
            
            eucl_dist  = nmf_euclidean_dist(X,W*H);
            eucl_dist_old  = nmf_euclidean_dist(X,Xr_old);
            errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
            fprintf('Iter %d ||X-WH||=%g, mean(X-WH)=%g, diff(k,k-1)=%g, ||X-WH||delta=%g\n',...
                iter,eucl_dist,errorx, diff, eucl_dist_old - eucl_dist);
            if errorx < 10^(-5), break, end
            
            Xr_old_models{model_iter} = Xr;
        end

         if record_scores
            diff_record(iter) = nmf_euclidean_dist(X,W*H);
            time_record(iter) = toc;
         end
         diff_record_models{model_iter} = diff_record;
    end
end

if (iter==maxiter)
    fprintf('max limit iteration reached (diff %g)\n', nmf_euclidean_dist(X,W*H) );
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
        otherwise 
            error('unknown solver - %s', als_solver);
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
