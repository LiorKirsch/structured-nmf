function  [best_W, best_H, best_diff_record, best_time_record, ...
          eucl_dist] = load_nmf_results(X, K, nmf_method, parms)
%
%
%

    structre_type = take_from_struct(parms, 'structre_type', 'none');
    
    vars = {'best_W', 'best_H', 'best_diff_record', ...
            'best_time_record','eucl_dist'};
    filename  = set_filenames('demixing', parms);
%     disp(filename);
    [do_calc, best_W, best_H, best_diff_record, best_time_record, ...
     eucl_dist] = cond_load(filename, 0, vars{1:end});
    if do_calc < 1 
       return
    end
    
    switch structre_type
        case 'tree', 
          tree_structure = parms.structure_matrix;
          tree_regions = parms.tree_regions;
          data_regions = parms.regions;
          [tree_X,reverse_map] = map_structure_to_region(X, data_regions, tree_regions);
          [tree_best_W, tree_best_H] = tree_nmf(tree_X, K, nmf_method, ...
              tree_structure,parms);
          best_H = tree_best_H(reverse_map);
          best_W = tree_best_W(reverse_map);
        case {'relations','relations_dist','relations_parentdist', 'relations_parent_level'}
%           X_models = split_data(X, sample_id);
          
          [best_W, best_H, best_diff_record, ...
              best_time_record, eucl_dist] ...
                    = nmf(X, K, 'alsWithRelations', parms);  
        case 'none', 
          [best_W, best_H, best_diff_record, ...
              best_time_record, eucl_dist] ...
                    = nmf(X, K, nmf_method, parms);  
        otherwise
            error('structure type is not supported: %s\n', parms.structre_type );
    end
        
    save(filename, vars{1:end});
    fprintf('Saved demixing results into [%s]\n', filename);
end

function X_array = split_data(X, sample_id)
    assert(size(X,1) == sample_id ,'each sample in X should have an id');
    uniq_samp_id = unique(sample_id);
    X_array = cell(length(uniq_samp_id) ,1);
    for i = 1:length(uniq_samp_id)
       X_array{i} = X( ismember( sample_id, uniq_samp_id(i) ) ,:);
    end

end

function [tree_X,reverse_map] = map_structure_to_region(X, data_regions, tree_regions)

    tree_X = cell(length(tree_regions),1);
    reverse_map = nan(length(data_regions),1);
    for tree_node_idx = 1:length(tree_regions)
        tree_region = tree_regions{tree_node_idx};
        
        data_regions_idx = find(ismember(data_regions,tree_region));
        
        if isempty(data_regions_idx)
            tree_X{tree_node_idx} = [];
        else
            tree_X{tree_node_idx} = X{data_regions_idx};
            reverse_map(data_regions_idx) = tree_node_idx;
        end
    end
    
end