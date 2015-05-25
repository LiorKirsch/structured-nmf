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
%    results.scores - mean correlation per hyperparameter combination
%    results.regions - the correpsonding regions  

  if length(loop_over_var_value) ==1
     % The stopping phase of the recursion
     region_names = parms.relation_regions;
     
     var_name = loop_over_var_name{1};
     var_values = loop_over_var_value{1};
     results = cell(length(var_values),1);
     
     % load the true profiles for each region
     mouse_cell_types = load('mouse_cell_type_profiles.mat');
     mouse_cell_types.expression = 2 .^ mouse_cell_types.expression;
     [true_celltype_per_region, celltypes_used] = match_region_with_true_profile(mouse_cell_types, region_names);
     
     % map the genes used in the true profile in the genes in the data
     [gene_inds_true_type,gene_inds_predictions] = compare_to_true_profile(...
         mouse_cell_types,gene_info, parms.species,parms);
      true_profiles = cellfun(@(x) x(gene_inds_true_type,:), ...
                         true_celltype_per_region,'UniformOutput',false);
                     
     parfor i_vars = 1: length(var_values)
         current_parms = parms;
         
         [new_loop_string,current_parms] = ...
             add_var_to_loop_string(loop_string, var_values, i_vars, ...
                                    var_name, current_parms);
         set_terminal_title(new_loop_string);
         
         [W, H] = load_nmf_results(X, current_parms.num_types, ...
                                   current_parms.nmf_method, current_parms);
         
         predicted_profiles = cellfun(@(x) x(:,gene_inds_predictions)', H,'UniformOutput',false);

         results{i_vars} = get_all_scores(predicted_profiles, W, ...
                                                        true_profiles, ...
                                                        region_names);  
         results{i_vars}.cell_type_used = celltypes_used;
         
     end
     
     baseline_celltype_profile = cellfun(@(x) mean(x), X,'UniformOutput',false);
     baseline_celltype_profile = cellfun(@(x) x(:,gene_inds_predictions)', baseline_celltype_profile,'UniformOutput',false);
     baseline_proportions = cellfun(@(x) zeros(3,3), X ,'UniformOutput',false);
     fprintf('\n===AVERAGE PROFILE SCORES===\n');
     get_all_scores(baseline_celltype_profile, baseline_proportions, ...
                                        true_profiles, region_names,false);  

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


function [new_loop_string,current_parms] = ...
        add_var_to_loop_string(loop_string, var_values, i_vars, ...
                               var_name, current_parms)
    if iscellstr(var_values)
       current_parms.(var_name) = var_values{i_vars};
       new_loop_string = sprintf('%s - %s %s',loop_string, var_name, ...
                                 current_parms.(var_name) );
   else
       current_parms.(var_name) = var_values(i_vars);
       new_loop_string = sprintf('%s - %s %g',loop_string, var_name, ...
                                 current_parms.(var_name) );
   end

end


function output = get_all_scores(predicted_profiles, predicted_proportions, ...
                                               true_profiles, ...
                                               region_name,show_indv)
   if ~exist('show_indv','var')
       show_indv = true;
   end
   
    num_regions = length(predicted_profiles);
    assert(length(true_profiles) == num_regions, ...
           'true and predicted profiles should have the same number of elements');        
    best_scores = nan(num_regions,1);
    indv_scores = cell(num_regions,1);
    for i = 1:num_regions
        % need to transpose the expression
        GT_proportions = zeros(size(predicted_proportions{i},1), ...
                               size(true_profiles{i},2));
        [~, ~, best_scores(i), ~,indv_scores{i}] = ...
            match_profiles_to_gt(predicted_proportions{i}, ...
                                 predicted_profiles{i}', true_profiles{i}', ...
                                 GT_proportions', 'spearman'); 
        fprintf('%s - %g\n',region_name{i}, best_scores(i));
        
        if show_indv
            for indv_i = 1:length(indv_scores{i})
                fprintf('\t\ttype %d - %g\n', indv_i, indv_scores{i}(indv_i));
            end
        end
    end
    output.scores =  best_scores;
    output.regions =  region_name;    
    
    output.mean_score = mean_corr_coeff(best_scores);
    fprintf('====REGIONS MEAN SCORE: %g====\n', output.mean_score);
end