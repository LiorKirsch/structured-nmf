function results = loopOverResults_real(X, gene_info,...
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
     
     var_name = loop_over_var_name{1};
     var_values = loop_over_var_value{1};
     results = cell(length(var_values),1);
     
     % for each region get only a subset of the genes
     [curr_gene_info, ~, parms] = gene_subset_selection(gene_info, ...
                                          X{1}, parms);
                                    
     for i_vars = 1: length(var_values)
         current_parms = parms;
         
         [new_loop_string, current_parms, var_current_value] = ...
             add_var_to_loop_string(loop_string, var_values, i_vars, ...
                                    var_name, current_parms);

         current_parms = get_current_markers(current_parms, curr_gene_info.gene_symbols);
         
      
         filename  = set_filenames('results', current_parms);
         load(filename, 'curr_result');
         results{i_vars} = curr_result;
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
          results{i_vars} = loopOverResults_real(X, gene_info, ...
                                                    current_parms, ...
                                                    loop_over_var_name, ...
                                                    loop_over_var_value, ...
                                                    new_loop_string);
      end 
  end
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
