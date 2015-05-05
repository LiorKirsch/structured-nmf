function  [best_W, best_H, best_diff_record, best_time_record, ...
          eucl_dist] = load_nmf_results(X, K, nmf_method, parms)
%
%
%

    if ~isfield(parms, 'structre_type')
        parms.structre_type = 'none';
    end
    
    vars = {'best_W', 'best_H', 'best_diff_record', ...
            'best_time_record','eucl_dist'};
    filename  = set_filenames('demixing', parms);

    [do_calc, best_W, best_H, best_diff_record, best_time_record, ...
     eucl_dist] = cond_load(filename, 0, vars{1:end});
    if do_calc < 1 
       return
    end
    
    switch parms.structre_type
        case 'tree', 
          tree_structure = parms.structure_matrix;
          tree_regions = parms.tree_regions;
          data_regions = parms.regions;
          [tree_X,reverse_map] = map_structure_to_region(X, data_regions, tree_regions);
          [tree_best_W, tree_best_H] = tree_nmf(tree_X, K, nmf_method, ...
              tree_structure,parms);
          best_H = tree_best_H(reverse_map);
          best_W = tree_best_W(reverse_map);
        case 'relations', 
          X_models = split_data(X, sample_id);
          relation_matrix_for_H = parms.structure_matrix;
          
          [best_W, best_H, best_diff_record, ...
              best_time_record, eucl_dist] ...
                    = nmf_als_with_relations(parms,X_models,...
                    relation_matrix_for_H,W_init_model, H_init_models);
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
    for i = 1:length(tree_regions)
        tree_region = tree_regions{i};
        
        z = find(ismember(data_regions,tree_region));
        
        if isempty(z)
            tree_X{i} = [];
        else
            tree_X{i} = X{z};
            reverse_map(z) = i;
        end
    end
    
end