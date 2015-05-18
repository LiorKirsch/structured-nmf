function [best_W,best_H,best_diff_record,best_time_record,eucl_dist] = nmf(X,K,alg,parms)
%
% NMF wrapper function
% function [W,H] = nmf(X,K,alg[,maxiter,speak])
%
% INPUT:
%           'X'     Inputmatrix - [num_samples X num_features]
%           'K'     Number of components
%           'alg'   Algorithm to use: 
%                   'mm'     multiplicative updates using euclidean
%                            distance. Lee, D..D., and Seung, H.S., (2001)
%                   'cjlin'  alternative non-negative least squares using 
%                            projected gradients, author: Chih-Jen Lin, 
%                            National Taiwan University.
%                   'prob'   probabilistic NFM interpretating X as samples
%                            from a multinomial, author: Lars Kai Hansen,
%                            Technical University of Denmark
%                   'als'    Alternating Least Squares. Set negative
%                            elements to zero. 
%                   'alsobs' Alternating Least Squares. Set negative elements
%                            to zero and adjusts the other elements acording
%                            to Optimal Brain Surgeon. 
%           'maxiter'   Maximum number of iterations, default = 1000.
%           'speak'     Print information to screen unless speak = 0,
%                       default = 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Based on the NMF-toolbox by Kasper Winther Joergensen
% http://cogsys.imm.dtu.dk/toolbox/nmf/
% 
% Lior Kirsch 03/2015


loglevel = take_from_struct(parms, 'loglevel', 1);
num_restarts = take_from_struct(parms, 'num_restarts', 1);
rand_seed = take_from_struct(parms, 'rand_seed', 42);
maxiter = take_from_struct(parms, 'maxiter', 1000);
do_sep_init = take_from_struct(parms, 'do_sep_init', false);

% find dimensionallity of X
[D,N] = size(X);
parms.debug  =1;
rng(rand_seed);
diff_record = nan(1,maxiter);
time_record = nan(1,maxiter);
eucl_dist = nan(num_restarts,1);
for i = 1:num_restarts
    
    if iscell(X) 
        [W_init, H_init] = cellfun(@(x) get_random_W_H(x,K,parms),X,'UniformOutput',false);
    else
        [W_init, H_init] = get_random_W_H(X,K,parms);
    end
    
    % switch algorithm 
    switch alg
        case 'mm'
            if loglevel, disp('Using mm algorithm'),end
            [W,H,diff_record,time_record]=nmf_mm(parms,X,W_init,H_init);
        case 'prob' 
            if loglevel, disp('Using prob algorithm'),end
            [W,H,diff_record,time_record]=nmf_prob(parms,X,W_init,H_init);
        case 'cjlin'
            if loglevel, disp('Using cjlin algorithm'),end
            [W,H,diff_record,time_record]=nmf_cjlin(parms,X,W_init,H_init);
        case 'alsPinv'
            parms.als_solver= 'pinv_project';
            if loglevel, disp('Using als pinv and project'),end
            [W,H,diff_record,time_record]=nmf_als(parms,X,W_init,H_init);
        case 'alsBlockpivot'
            parms.als_solver= 'blockpivot';
            if loglevel, disp('Using als-blockpivot algorithm'),end
            [W,H,diff_record,time_record]=nmf_als(parms,X,W_init,H_init);
        case 'alsActiveSet'
            parms.als_solver= 'active_set';
            if loglevel, disp('Using als-active-set algorithm'),end
            [W,H,diff_record,time_record]=nmf_als(parms,X,W_init,H_init);
        case 'alsWithRelations'
            switch parms.nmf_method
                case 'alsActiveSet'
                    parms.als_solver= 'active_set';
                case 'alsBlockpivot'
                    parms.als_solver= 'blockpivot';
                case 'alsPinv'
                    parms.als_solver= 'pinv_project';
                otherwise
                    error('Unknown nmf method %s', parms.nmf_method);
            end
            
            if do_sep_init
                curr_parms = parms;
                curr_parms.H_lambda = 0;
                [W_init,H_init]=cellfun(@(x,w,h) nmf_als(curr_parms,x,w,h) ,X,W_init,H_init,'UniformOutput',false);
            end
            
            if isinf(parms.H_lambda)
                % use samples from all X{i}, build a tree that is 
                curr_parms = parms;
                curr_parms.H_lambda  = 0;
                curr_parms.W_lambda  = 0;
                
                reverse_map = [];
                
                dim = ndims(H_init{1});          % Get the number of dimensions for your arrays
                M = cat(dim+1,H_init{:});        % Convert to a (dim+1)-dimensional matrix
                curr_H_init = mean(M,dim+1);     % Get the mean across arrays
                
                curr_W_init = cat(1,W_init{:});  % Concat W_init
                curr_X = cat(1,X{:});            % Concat X
                for i_nodes = 1:length(X)
%                     curr_X = cat(1,curr_X,X{i_nodes});
                    reverse_map = cat(1,reverse_map, i_nodes*ones(size(X{i_nodes},1),1) );
                end
                
                [W_combined,H_combined] = nmf_als(curr_parms,curr_X,curr_W_init,curr_H_init);
                for i_nodes =1:length(X)
                    H{i_nodes} = H_combined;
                    W{i_nodes} = W_combined(reverse_map==i_nodes,:);
                end
            else
                relation_matrix_for_H = parms.structure_matrix;
                if loglevel, disp('Using als-with-relations'),end
                [W,H,diff_record,time_record]=nmf_als_with_relations(parms,X,...
                        relation_matrix_for_H,W_init, H_init);
            end
    %     case 'alsobs'
    %         if loglevel, disp('Using alsobs algorithm'),end
    %         [W,H]=nmf_alsobs(X,K,maxiter,loglevel);
        otherwise
            error('Unknown method. Type "help nmf" for usage.');
            return
    end

    if iscell(X) && iscell(W) && iscell(H)
       num_elements = length(X);
       assert(length(W) == num_elements, ' X and W should have the same number of elements');
       assert(length(H) == num_elements, ' X and H should have the same number of elements');

       err = nan(num_elements,1);
       for m=1:num_elements
          err(m) =  nmf_euclidean_dist(X{m},W{m}*H{m});
       end 
       eucl_dist(i) = sum(err); % compute the sum of err over components
    else
       eucl_dist(i) = nmf_euclidean_dist(X,W*H);
    end
    
    if i==1
        best_W = W;
        best_H = H;
        best_diff_record = diff_record;
        best_time_record = time_record;
        best_iteration = 1;
    else
        if eucl_dist(i) < eucl_dist(i-1)
            best_H = H;
            best_W = W;
            best_diff_record = diff_record;
            best_time_record = time_record;
            best_iteration = i;
        end
    end
end

fprintf('==== best random restart - %d ====\n',best_iteration);
end
