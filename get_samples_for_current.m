function [curr_X,reverse_map] = get_samples_for_current(X,current_id, node_to_include,parms)
    
    curr_X = [];
    reverse_map = [];
    fprintf('region %s contains:\n', parms.tree_regions{current_id} );
    for i = 1:length(node_to_include)
        child_id = node_to_include(i);
        curr_X = cat(1,curr_X,X{child_id});
        reverse_map = cat(1,reverse_map, child_id*ones(size(X{child_id},1),1) );
        fprintf('\t%s (%d)\n', parms.tree_regions{child_id}, size(X{child_id},1) );
    end
    
end