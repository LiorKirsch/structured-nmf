function [W_nodes,H_nodes] = tree_nmf(X ,K ,alg , tree_structure,parms)
% When you have a hierarcy of regions
% first learn a factorization that is good for the coarse region.
% then for each child region learn a unique factorization that uses the
% coarse region as prior.
%
% X - matrix of size (num_sample X num_features)
% K - number of components in the factorization
% alg - the nmf solver
% sample_node_id - for each sample specify the node in the tree [1,num_nodes]
% tree_structure - [num_nodes X num_nodes] holds the parent-child reletions
%    where tree_structure(i,j)=1 means i is a parent of j
% parms - hyperparms
%   parms.H_lambda - controls the strength of the H connection in the tree
%    (default 0 - connection does not effect the strength)
%   parms.W_lambda - controls the strength of the W connection in the tree
%    (default 0 - connection does not effect the strength)
%
% TODO:Baseline to compare to
%      a model that is unique for each child node
%      a model that is shared across all regions
%      a model that share the profiles but solves a least squares
%        problem for each child without taking the relation into account
%      
%      the 0 and inf case only support H tree (not W tree)
% Current version supports at most one parent per node (tree structure).
%
assert( parms.W_lambda == 0,'only H structure is supported so far');
    num_nodes = size(tree_structure,1);
assert( length(X) == num_nodes,'every node in the tree should be have a samples in X');

    H_nodes = cell(num_nodes,1);
    W_nodes = cell(num_nodes,1);

    
    if (parms.H_lambda == inf) 
        % use samples from all X{i}, build a tree that is 
        parms.H_lambda  = 0;
        parms.W_lambda  = 0;
        [curr_X,reverse_map] = get_samples_for_current(X,1, (1:num_nodes),parms);
        [W,H,~,~,~] = nmf(curr_X,K,alg,parms);
        for i =1:num_nodes
            H_nodes{i} = H;
            W_nodes{i} = W(reverse_map==i,:);
        end
    else   
       if (parms.H_lambda == 0) % only use leafs - don't use the tree structure
           tree_structure = false(num_nodes);
       end
       

       child_level_node_inds = find(sum(tree_structure,2) ==0) ;

       for tree_i = 1:length(child_level_node_inds)
            current_node_id = child_level_node_inds(tree_i);
            [H,W, H_nodes,W_nodes] = recursive_do_nmf(X ,K ,alg , tree_structure,current_node_id, H_nodes, W_nodes,parms);
            H_nodes{current_node_id} = H;
            W_nodes{current_node_id} = W;
       end
    end
    
    
end

function [H,W, H_nodes,W_nodes] = recursive_do_nmf(X ,K ,alg , tree_structure, structure_node_id, H_nodes, W_nodes, parms)

    node_parents = find( tree_structure(:,structure_node_id));
    
    
    % recursion step: do nmf for each parent
    for i = 1:length(node_parents)
       parent_node_id = node_parents(i);
       if  isempty( H_nodes{ parent_node_id } )
           % calc nmf for parent
           [H_parent,W_parent, H_nodes,W_nodes] = recursive_do_nmf(X ,K ,alg , tree_structure, parent_node_id, H_nodes, W_nodes, parms);
           H_nodes{parent_node_id} = H_parent;
           W_nodes{parent_node_id} = W_parent;
       end
    end
    
    
    %once you are done with all of the parents calc the current node
    % recursion stop (no parents OR all parents already computed)
    mean_parent_H = mean(cat(3,H_nodes{node_parents}),3);
    mean_parent_W = mean(cat(3,W_nodes{node_parents}),3);
    
    current_parms = parms;
    
    if isempty(node_parents)
        current_parms.H_lambda  = 0;
        current_parms.W_lambda  = 0;
    else
       %TODO USE the tree to decide a new H_lambda 
    end
    current_parms.H_prior = mean_parent_H;
    current_parms.W_prior = mean_parent_W;
    
    node_child_recursive =  inv(eye(size(tree_structure)) - tree_structure); % including self
    node_child_recursive = find( node_child_recursive(structure_node_id,:) );
    curr_X = get_samples_for_current(X,structure_node_id, node_child_recursive,parms);
    if isempty(curr_X)
        W = [];
        H = [];
    else
        [W,H,~,~,~] = nmf(curr_X,K,alg,current_parms);
    end
end


