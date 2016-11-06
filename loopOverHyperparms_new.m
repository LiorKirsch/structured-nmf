function [results, H_results, W_results,parms_results] = loopOverHyperparms_new(parms, loop_over_var_name, ...
                                     loop_over_var_value , loop_type)
% This function calls the nmf factorization but with different parms
% it loops over the set of parms provided in 
%  loop_over_var_name
%     and
%  loop_over_var_value
%
%
% The output is a strcuture that contains 
%    results.individual_scores =  score for each of different cell types
%           K scores for each region.
%    results.scores - the average score over the celltype (after z transform)
%           1 score for each region.
%    results.regions - the correpsonding regions  
%    results.mean_score - the average score (across all regions) - 1 score.
 

  [parms_list, parms_strings] = get_parms_list(parms, loop_over_var_name, loop_over_var_value);
  num_parms_comb = length(parms_list);
  
  results = cell(num_parms_comb,1);
  H_results = cell(num_parms_comb,1);
  W_results = cell(num_parms_comb,1);
  parms_results = cell(num_parms_comb,1);
  
  if strcmp(loop_type, 'rand_run')
      fprintf('shuffling the random seed\n');
      rng('shuffle');
      rand_perm = randperm( num_parms_comb );
      parms_list = parms_list( rand_perm );
      parms_strings = parms_strings( rand_perm );
  end
  
  for i =1 : num_parms_comb
    
    current_parms = parms_list{i};
    new_loop_string = parms_strings{i};
    set_terminal_title(new_loop_string);
    
    if current_parms.is_synthetic
        mix_data = create_multi_region_mix(fullfile(current_parms.mix_dir, current_parms.mix_files));

        mix_data.profiles = mix_data.profiles{1};
        mix_data.proportions = mix_data.proportions{1};
        mix_data.expression = mix_data.expression{1};
        mix_data.cell_types = mix_data.cell_types{1};
        mix_data.region = mix_data.region{1};
        
        mix_data.expression = change_to_linear_scale(mix_data.expression);
        
        curr_X = mix_data.expression;
        true_profiles = mix_data.profiles;
        
        current_parms.regions = mix_data.region;
        current_parms.cell_types = mix_data.cell_types;
    else
        [current_parms.dataset_file, expression, gross_region_vec, gene_info, ...
            gross_structures_info] = load_all(current_parms.dataset, current_parms.regions_ack);
        expression = change_to_linear_scale(expression);
        gross_regions = gross_structures_info(gross_region_vec);
        [X, region_names] = split_to_cell(expression, gross_regions);

        current_parms = get_structure_matrix(current_parms.dataset_file, current_parms.structre_type,region_names, current_parms);
        
        [current_parms.gene_hash, region_names, curr_gene_info,curr_X, gene_inds_true_type,...
          gene_inds_predictions, true_profiles, baseline_result ,...
          randbase_result, mouse_cell_types, celltypes_used] = ...
          load_cached_baseline_and_genesubsets(...
          current_parms, X, gene_info);
    end
          
    switch loop_type
        case {'rand_run', 'ordered_run'}

            load_nmf_results(curr_X, current_parms.num_types, ...
                                   current_parms.nmf_method, current_parms);
          
        case 'collect_results'

    %          current_parms = get_current_markers(current_parms, curr_gene_info.gene_symbols);

             [W, H] = load_nmf_results(curr_X, current_parms.num_types, ...
                                       current_parms.nmf_method, current_parms);

             predicted_profiles = cellfun(@(x) x(:,gene_inds_predictions)', ...
                                          H, 'UniformOutput', false);

             curr_result = get_all_scores(predicted_profiles, W, ...
                                                            true_profiles, ...
                                                            region_names);
             curr_result.parms = current_parms;
             curr_result.cell_type_used = celltypes_used;
             curr_result.loop_string = new_loop_string;
             
             curr_result.baseline_score = baseline_result.run_score;
             curr_result.baseline_celltype_region_avg_scores = baseline_result.celltype_region_avg_scores;
             curr_result.baseline_region_scores = baseline_result.region_scores ;
             curr_result.baseline_celltype_score = baseline_result.celltype_scores ;

             curr_result.randbase_score = randbase_result.run_score;
             curr_result.randbase_score_sem = randbase_result.run_score_sem;
             curr_result.randbase_celltype_region_avg_scores = randbase_result.celltype_region_avg_scores;
             curr_result.randbase_celltype_region_avg_scores_sem = randbase_result.celltype_region_avg_scores_sem;
             curr_result.randbase_region_scores = randbase_result.region_scores ;
             curr_result.randbase_region_scores_sem = randbase_result.region_scores_sem ;
             curr_result.randbase_celltype_score = randbase_result.celltype_scores ;

             parsave(curr_result,current_parms);
             results{i} = curr_result;
             H_results{i} = H;
             W_results{i} = W;
             parms_results{i} = current_parms;
         
        otherwise
            erorr('unkown loop_type')
    end
  end
  
  set_terminal_title('done');

end



function output = get_all_scores(predicted_profiles, predicted_proportions, ...
                                               true_profiles, ...
                                               region_name, show_individual)
   if ~exist('show_individual','var')
       show_individual = true;
   end
   
    num_regions = length(predicted_profiles);
    assert(length(true_profiles) == num_regions, ...
           'true and predicted profiles should have the same number of elements');        
    best_scores = nan(num_regions,1);
    individual_scores = cell(num_regions,1);

    for i = 1:num_regions
        % need to transpose the expression
        GT_proportions = zeros(size(predicted_proportions{i},1), ...
                               size(true_profiles{i},2));
        [~, ~, best_scores(i), ~,individual_scores{i}] = ...
            match_profiles_to_gt(predicted_proportions{i}, ...
                                 predicted_profiles{i}', true_profiles{i}', ...
                                 GT_proportions', 'spearman'); 
        
        % when the number of true cell type is less then 3 pad the vector with nans
        individual_scores{i} = [individual_scores{i}; ...
                nan( max(0,3 -length(individual_scores{i})),1)]; 
            
        %         fprintf('%25s: %4.2f ',region_name{i}, 100*best_scores(i));%         
        %         if show_individual
        %             for individual_i = 1:length(individual_scores{i})
        %                 fprintf('\ttype %d: %4.2f', individual_i, 100*individual_scores{i}(individual_i));
        %             end
        %             fprintf('\n');
        %         end
    end
    output.celltype_scores =  individual_scores;
    output.region_scores =  best_scores;
    output.regions =  region_name;        
    output.run_score = mean_corr_coeff(best_scores);

    dim = ndims(individual_scores);          % Get the number of dimensions for your arrays
    M = cat(dim+1,individual_scores{:});        % Convert to a (dim+1)-dimensional matrix
    output.celltype_region_avg_scores = nanmean(M,dim+1);     % Get the mean across arrays
                
    %     fprintf('\tREGION MEAN SCORE: %5.3g====\n', output.run_score);
end


function parsave(curr_result,parms)
%
    filename  = set_filenames('results', parms);
    save(filename, 'curr_result','parms');
end
