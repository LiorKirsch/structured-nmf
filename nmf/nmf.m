function [best_W,best_H,best_diff_record,best_time_record,eucl_dist] ...
        = nmf(X, K, alg, parms)
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


   % loglevel = take_from_struct(parms, 'loglevel', 1);
   num_restarts = take_from_struct(parms, 'num_restarts', 1);
   init_type = take_from_struct(parms, 'init_type', 'random');
   % maxiter = take_from_struct(parms, 'maxiter', 1000);
   % do_sep_init = take_from_struct(parms, 'do_sep_init', false);
   
   % find dimensionallity of X
   % [D, N] = size(X);
   parms.debug  =1;
 
   H = cell(1, num_restarts);   
   diff_record = H; time_record = H; W = H;
   eucl_dist = nan(num_restarts,1);
   
   % finding genes with zero expression
   
   
   if iscell(X) 
       num_genes = size(X{1},2);
       genes_with_zero_expression = true(1,num_genes);
       for j_regions = 1:length(X);
           genes_with_zero_expression = genes_with_zero_expression & ...
               all(X{j_regions} ==0,1) ;
       end
       % I remove gene which have zero expression no and return them at the end
       X = cellfun(@(x) x(:,~genes_with_zero_expression), X,'UniformOutput',false);
   else
      num_genes = size(X,2);
      genes_with_zero_expression = all(X==0,1) ;
      X = X(:,~genes_with_zero_expression); 
   end
   
   switch init_type
     case 'random'
       parfor i = 1:num_restarts
           parms_current = parms;
           parms_current.restart_ind = i;
           parms_current = rmfield(parms_current,'num_restarts'); % so we can share the restarts even when the total restart changes
           
           fprintf('====== restart %d (%d) HL = %4.2g\n', i, ...
                   num_restarts, parms_current.H_lambda);
           
            filename  = set_filenames('demixing_rand_restart', parms_current);
           vars = {'current_W', 'current_H', 'current_diff_record', ...
            'current_time_record','current_eucl_dist'};
        
            [do_calc, current_W, current_H, current_diff_record, current_time_record, ...
             current_eucl_dist] = cond_load(filename, 0, vars{1:end});

            if do_calc < 1 
               fprintf('loading random restart from memory - %s\n', filename);
               % do thing
            else
               if iscell(X) 
                   [W_init, H_init] = cellfun(@(x) get_random_W_H(x, K, ...
                                                                     parms_current), X, 'UniformOutput', false);
               else
                   [W_init, H_init] = get_random_W_H(X,K,parms_current);
               end

               [current_W, current_H, current_diff_record, current_time_record] = ...
                   nmf_alg_selection(X,W_init,H_init,alg,parms_current);           
               current_eucl_dist = compute_eucl_dist(X, current_W,current_H);
               
               
               parsave(filename, current_W, current_H, current_diff_record,...
                   current_time_record, current_eucl_dist);
               fprintf('Saved demixing random restart model into [%s]\n', filename);
            end
            W{i} = current_W;
            H{i} = current_H;
            diff_record{i} = current_diff_record;
            time_record{i} = current_time_record;
            eucl_dist(i) = current_eucl_dist;
            
       end

       % Select the best 
       [~, best_iteration] = min(eucl_dist);
       best_H = H{best_iteration};
       best_W = W{best_iteration};
       best_diff_record = diff_record{best_iteration};
       best_time_record = time_record{best_iteration};
       fprintf('=== nmf: best restart = %d  HL=%4.2g , best-worst=%4.2g\n', best_iteration, ...
               parms.H_lambda , max(eucl_dist) - min(eucl_dist) );
       
     case 'svd'
       if iscell(X) 
           [W_init, H_init] = cellfun(@(x) nndsvd(x,K,0),X,'UniformOutput',false);
       else
           [W_init, H_init] = nndsvd(X, K,0);
       end
       [best_W,best_H,best_diff_record,best_time_record] = nmf_alg_selection(X,W_init,H_init,alg,parms);
       eucl_dist = compute_eucl_dist(X,best_W,best_H);
     otherwise
       error('unkown init option -%s',init_type);    
   end
   
   
   % I now return genes with zero expression 
   if iscell(X) 
       for i = 1:length(X)
            new_H = zeros(K, num_genes);
            new_H(:,~genes_with_zero_expression) = best_H{i};
            best_H{i} = new_H;
       end
   else
       new_H = zeros(K, num_genes);
       new_H(:,~genes_with_zero_expression) = best_H;
       best_H = new_H;
   end
end


function parsave(filename,current_W, current_H, current_diff_record,...
                   current_time_record, current_eucl_dist)
%
    save(filename, 'current_W', 'current_H', 'current_diff_record',...
                   'current_time_record', 'current_eucl_dist');
end

function parsave_warm(filename,W_warm, H_warm)
%
    save(filename, 'W_warm', 'H_warm');
end

function [W, H, diff_record, time_record] = nmf_alg_selection(X,W_init,H_init,alg,parms)
    loglevel = take_from_struct(parms, 'loglevel', 1);
    do_sep_init = take_from_struct(parms, 'do_sep_init', false);
    diff_record = nan;
    time_record = nan;

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

            %%% TODO hot start with good HL=0.
             if do_sep_init
%                 fprintf('=== warm start === (HL=%4.2g)\n', parms.H_lambda);
                warm_parms = parms;
                warm_parms.H_lambda = 0; % !!!! must set lambda to zero 0!!!!
                
                
                 [~, warm_filename, warm_dir]  = set_filenames('demixing_rand_restart', warm_parms);
                 warm_filename = fullfile(warm_dir,['warm', warm_filename]);
                   vars = {'W_warm', 'H_warm'};
        
                [do_calc, W_warm, H_warm ]= cond_load(warm_filename, 0, vars{1:end});

                if do_calc < 1 
                   fprintf('loading warm restart from memory %s\n', warm_filename);
                else
                    for i_cell =1 :length(X)
                        if isfield(warm_parms,'H_markers')
                            warm_parms.H_markers = parms.H_markers{i_cell};
                        end
                        [W_warm{i_cell}, H_warm{i_cell}, ~, ~] = ...
                            nmf_als(warm_parms, X{i_cell}, ...
                                                W_init{i_cell}, H_init{i_cell});
                    end
                    fprintf('Saving warm restart %s\n', warm_filename);
                    parsave_warm(warm_filename,W_warm, H_warm);
                    
                end
                W_init = W_warm;
                H_init = H_warm;
                
                % [W_init,H_init]=cellfun(@(x,w,h) nmf_als(curr_parms,x,w,h) ,X,W_init,H_init,'UniformOutput',false);
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
                
                if isfield(curr_parms,'H_markers')
                    dim = ndims(curr_parms.H_markers{1}); % Get the number of dimensions for your arrays
                    M = cat(dim+1,curr_parms.H_markers{:}); % Convert to a (dim+1)-dimensional matrix
                    curr_parms.H_markers = any(M,dim+1); % Get any marker
                end
                
                curr_W_init = cat(1,W_init{:});  % Concat W_init
                curr_X = cat(1,X{:});            % Concat X
                for i_nodes = 1:length(X)
                    % curr_X = cat(1,curr_X,X{i_nodes});
                    reverse_map = cat(1,reverse_map, i_nodes*ones(size(X{i_nodes},1),1) );
                end
                
                [W_combined, H_combined] = nmf_als(curr_parms, ...
                                                   curr_X, ...
                                                   curr_W_init,curr_H_init);
                H = cell(1, length(X));
                W = cell(1, length(X));
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
            % case 'alsobs'
            %   if loglevel, disp('Using alsobs algorithm'),end
            %   [W,H]=nmf_alsobs(X,K,maxiter,loglevel);

      otherwise
            error('Unknown method. Type "help nmf" for usage.');
    end

end


function eucl_dist = compute_eucl_dist(X, W, H)

    if iscell(X) && iscell(W) && iscell(H)
       num_elements = length(X);
       assert(length(W) == num_elements, ' X and W should have the same number of elements');
       assert(length(H) == num_elements, ' X and H should have the same number of elements');

       err = nan(num_elements,1);
       for m=1:num_elements
          err(m) =  nmf_euclidean_dist(X{m},W{m}*H{m});
       end 
       eucl_dist = sum(err); % compute the sum of err over components
    else
       eucl_dist = nmf_euclidean_dist(X,W*H);
    end
                
end