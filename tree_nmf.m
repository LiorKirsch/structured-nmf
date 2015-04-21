function tree_nmf(X ,K ,alg , sample_node_id, tree_structure,parms)
% When you have a hierarcy of regions
% first learn a factorization that is good for the coarse region.
% then for each child region learn a unique factorization that uses the
% coarse region as prior.
%
% X - matrix of size (num_sample X num_features)
% K - number of components in the factorization
% alg - the nmf solver
% sample_node_id - for each sample specify the node in the tree [1,num_nodes]
% tree_structure - [num_nodes X num_nodes] matrix which hold father child
%     reletions where tree_structure(i,j)=1 means i is a parent of j
% parms - hyperparms
%   parms.H_lambda - controls the strength of the H connection in the tree
%    (default 0 - connection does not effect the strength)
%   parms.W_lambda - controls the strength of the W connection in the tree
%    (default 0 - connection does not effect the strength)
%
% TODO:Baseline to compare to
%      a model that is unique for each child node
%      a model that is shared across all region
%      a model that share the profiles but than solve a least squares
%        problem for each child
% 
% Current version support at most one parent per node.

    num_nodes = size(tree_structure,1);
    H_nodes = cell(num_nodes,1);
    W_nodes = cell(num_nodes,1);

    child_level_node_inds = find(sum(tree_structure,2) ==0) ;

    for i = 1:length(child_level_node_inds)
        current_node_id = child_level_node_inds(i);
        [H,W, H_nodes,W_nodes] = recursive_do_nmf(current_node_id, tree_structure, H_nodes, W_nodes);
        H_nodes{i} = H;
        W_nodes{i} = W;
    end

end

function [H,W, H_nodes,W_nodes] = recursive_do_nmf(node_id, tree_structure, H_nodes, W_nodes)

    node_parents = find( tree_structure(:,node_id));
    
    
    % recursion step: do nmf for each parent
    for i = 1:node_parents
       parent_node_id = node_parents(i);
       if  isempty( H_nodes{ parent_node_id } )
           % calc nmf for parent
           [H_parent,W_parent] = recursive_do_nmf(parent_node_id, tree_structure, H_nodes, W_nodes);
           H_nodes{parent_node_id} = H_parent;
           W_nodes{parent_node_id} = W_parent;
       end
    end
    
    %once you are done with all of the parents calc the current node
    % recursion stop (no parents OR all parents already computed)
    mean_parent_H = mean(cat(3,H_nodes{node_parents}),3);
    mean_parent_W = mean(cat(3,W_nodes{node_parents}),3);
    
    current_parms = parms;
    current_parms.H_prior = mean_parent_H;
    current_parms.W_prior = mean_parent_W;
    
    [W,H,~,~,~] = nmf(X,K,alg,parms);
end