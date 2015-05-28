function results = loopOverHyperparms_real(X, gene_info,...
                                     parms, loop_over_var_name, ...
                                     loop_over_var_value ,loop_string)
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
    

  if length(loop_over_var_value) ==1
     % The stopping phase of the recursion
     region_names = parms.relation_regions;
     
     var_name = loop_over_var_name{1};
     var_values = loop_over_var_value{1};
     results = cell(length(var_values),1);
     
     % for each region get only a subset of the genes
     [curr_gene_info, ~, parms] = gene_subset_selection(gene_info, ...
                                          X{1}, parms);
     
     [~, curr_X, ~] = cellfun(@(x) gene_subset_selection(gene_info, ...
                                          x, parms), X,'UniformOutput',false);
     
   
    
     % load the true profiles for each region
     mouse_cell_types = load('mouse_cell_type_profiles.mat');
     mouse_cell_types.expression = 2 .^ mouse_cell_types.expression;
     [true_celltype_per_region, celltypes_used] = ...
         match_region_with_true_profile(mouse_cell_types, region_names);
     
     % map the genes used in the true profile in the genes in the data
     [gene_inds_true_type,gene_inds_predictions] = compare_to_true_profile(...
         mouse_cell_types,curr_gene_info, parms.species,parms);
      true_profiles = cellfun(@(x) x(gene_inds_true_type,:), ...
                         true_celltype_per_region,'UniformOutput',false);
         
     parfor i_vars = 1: length(var_values)
         current_parms = parms;
         
         [new_loop_string, current_parms, var_current_value] = ...
             add_var_to_loop_string(loop_string, var_values, i_vars, ...
                                    var_name, current_parms);
         set_terminal_title(new_loop_string);

         current_parms = get_current_markers(current_parms, curr_gene_info.gene_symbols);
         
         [W, H] = load_nmf_results(curr_X, current_parms.num_types, ...
                                   current_parms.nmf_method, current_parms);
         
         predicted_profiles = cellfun(@(x) x(:,gene_inds_predictions)', ...
                                      H, 'UniformOutput', false);

         results{i_vars} = get_all_scores(predicted_profiles, W, ...
                                                        true_profiles, ...
                                                        region_names);
         results{i_vars}.cell_type_used = celltypes_used;
         results{i_vars}.var_name = var_name;
         results{i_vars}.loop_string = new_loop_string;
         results{i_vars}.var_value = var_current_value;
         
     end
     
     baseline_celltype_profile = cellfun(@(x) mean(x), curr_X,'UniformOutput',false);
     baseline_celltype_profile = cellfun(@(x) x(:,gene_inds_predictions)', ...
                                         baseline_celltype_profile,'UniformOutput',false);
     baseline_celltype_profile = cellfun(@(x) repmat(x,1,3), ...
                                         baseline_celltype_profile,'UniformOutput',false);
     baseline_proportions = cellfun(@(x) zeros(3,3), curr_X ,'UniformOutput',false);
%      fprintf('\n===AVERAGE PROFILE SCORES===\n');
     baseline_result = get_all_scores(baseline_celltype_profile, baseline_proportions, ...
                                        true_profiles, region_names,false);  
     for i_vars = 1: length(var_values)
         results{i_vars}.baseline_score = baseline_result.run_score;
         results{i_vars}.baseline_celltype_score = baseline_result.celltype_region_avg_scores;
     end
  else
      % Remove one layer from the recursion
      var_name = loop_over_var_name{1};
      loop_over_var_name = loop_over_var_name(2:end);
      var_values = loop_over_var_value{1};
      loop_over_var_value = loop_over_var_value(2:end);
      
      results = cell(length(var_values),1);
      for i_vars = 1: length(var_values)
          current_parms = parms;
          [new_loop_string,current_parms] = ...
              add_var_to_loop_string(loop_string, var_values, i_vars, ...
                                     var_name, current_parms);
          results{i_vars} = loopOverHyperparms_real(X, gene_info, ...
                                                    current_parms, ...
                                                    loop_over_var_name, ...
                                                    loop_over_var_value, ...
                                                    new_loop_string);
      end 
  end
  set_terminal_title('done');
end


function [new_loop_string, current_parms, value] = ...
        add_var_to_loop_string(loop_string, var_values, i_vars, ...
                               var_name, current_parms)
    if iscellstr(var_values)
       value = var_values{i_vars};
       current_parms.(var_name) = value;
       new_loop_string = sprintf('%s - %s %s', loop_string, var_name, ...
                                 current_parms.(var_name) );
    else
       value = var_values(i_vars);
       current_parms.(var_name) = value;
       new_loop_string = sprintf('%s - %s %5.3g', loop_string, var_name, ...
                                 current_parms.(var_name) );
   end

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
%         fprintf('%25s: %4.2f ',region_name{i}, 100*best_scores(i));
%         
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
    output.celltype_region_avg_scores = mean(M,dim+1);     % Get the mean across arrays
                
%     fprintf('\tREGION MEAN SCORE: %5.3g====\n', output.run_score);
end

function parms = get_current_markers(parms, gene_symbols)

    num_genes = length(gene_symbols);
    num_regions = length(parms.relation_regions);
     if isfield(parms,'num_markers')
        [neuro_mrk,astro_mrk,oligo_mrk] = get_okaty_markers(parms.num_markers, 1000);
        neuro_inds = get_intersecting_genes(...
            gene_symbols, neuro_mrk, parms);
        astro_inds = get_intersecting_genes(...
            gene_symbols, astro_mrk, parms);
        oligo_inds = get_intersecting_genes(...
            gene_symbols, oligo_mrk, parms);
        
        H_markers = false(parms.num_types,num_genes);
        H_markers(1, : ) = neuro_inds;
        H_markers(2, :) = astro_inds;
        H_markers(3, : ) = oligo_inds;

        % TODO - change so each region has its own markers
        parms.H_markers = repmat({H_markers}, num_regions,1);
    end
    
end