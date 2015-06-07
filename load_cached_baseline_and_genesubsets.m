function [gene_hash, region_names, curr_gene_info,curr_X, gene_inds_true_type,...
          gene_inds_predictions, true_profiles, baseline_result ,...
          randbase_result, mouse_cell_types, celltypes_used] = ...
          load_cached_baseline_and_genesubsets(...
          parms, X, gene_info)

        cache_parms = parms;
        cache_parms = rmfield(cache_parms,'H_lambda');
        cache_parms = rmfield(cache_parms,'W_lambda');
        cache_parms = rmfield(cache_parms,'init_type');
        cache_parms = rmfield(cache_parms,'init_subtype');
        cache_parms = rmfield(cache_parms,'W_constraints');
        cache_parms = rmfield(cache_parms,'maxiter');
        cache_parms = rmfield(cache_parms,'do_sep_init');
        cache_parms = rmfield(cache_parms,'num_restarts');
        cache_parms = rmfield(cache_parms,'nmf_method');
        
        true_dataset = take_from_struct(cache_parms, 'true_dataset', 'okaty');   %'barres'  


        
       filename  = set_filenames('baselines', cache_parms);
           vars = {'gene_hash', 'curr_gene_info', 'curr_X', ...
                   'gene_inds_true_type','gene_inds_predictions'...
                   'true_profiles','baseline_result'...
                   'randbase_result','mouse_cell_types'...
                   'celltypes_used','region_names'};
        
            [do_calc, gene_hash, curr_gene_info, curr_X, gene_inds_true_type, ...
             gene_inds_predictions, true_profiles, baseline_result, randbase_result, ...
             mouse_cell_types, celltypes_used,region_names] = cond_load(filename, 0, vars{1:end});
         
        if do_calc < 1 
               fprintf('loading baselines from cache - %s\n', filename);
            else

         % for each region get only a subset of the genes
             region_names = parms.relation_regions;

             [curr_gene_info, ~, gene_hash] = gene_subset_selection(gene_info, ...
                                                               X{1}, cache_parms);

             [~, curr_X, ~] = cellfun(@(x) gene_subset_selection(gene_info, ...
                                                  x, cache_parms), X,'UniformOutput',false);

             % load the true profiles for each region
              switch true_dataset
                 case 'barres'
                    mouse_cell_types = load_data('barres2014');
                    mouse_cell_types.expression = mouse_cell_types.data;
                    mouse_cell_types.all_symbols = mouse_cell_types.gene_symbols;
                    mouse_cell_types.gene_symbol = mouse_cell_types.gene_symbols;
                    mouse_cell_types.refer_to_index = 1:length(mouse_cell_types.gene_symbol);
                  case 'okaty'
                      mouse_cell_types = load('mouse_cell_type_profiles.mat');
                      mouse_cell_types.expression = 2 .^ mouse_cell_types.expression;
                  otherwise
                      error('unkown ground true dataset - %s', true_dataset)
              end
             
             [true_celltype_per_region, celltypes_used] = ...
                 match_region_with_true_profile(mouse_cell_types, region_names, parms);

             % map the genes used in the true profile in the genes in the data
             [gene_inds_true_type,gene_inds_predictions] = compare_to_true_profile( ...
                 mouse_cell_types,curr_gene_info, cache_parms.species,cache_parms);
              true_profiles = cellfun(@(x) x(gene_inds_true_type,:), ...
                                      true_celltype_per_region,'UniformOutput',false);

             %===== mean profile baseline ======%
             [baseline_celltype_profile, baseline_proportions] = get_mean_baseline( ...
                 curr_X,gene_inds_predictions);
             fprintf('loopoverhyperparms_real: create mean baseline\n')
             baseline_result = get_all_scores(baseline_celltype_profile, ...
                                              baseline_proportions, ...
                                              true_profiles, region_names,false);  

             %===== rand sample baseline ======%
             %      [randbase_celltype_profile, randbase_proportions] =...
             %            get_randsamp_baseline(curr_X,gene_inds_predictions,parms);
             %      randbase_result = get_all_scores(randbase_celltype_profile, randbase_proportions, ...
             %                                         true_profiles, region_names,false);           

             fprintf('loopoverhyperparms_real: create random baseline\n');
             randbase_result = get_randbaseline(curr_X, ...
                                                gene_inds_predictions,region_names, ...
                                                true_profiles, ...
                                                cache_parms.num_types,cache_parms);
             save(filename, vars{1:end});
        end
        
end