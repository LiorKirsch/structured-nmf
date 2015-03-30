function scores = loopOverHyperparms(X,GT_profiles,rand_subset,parms, loop_over_var_name, loop_over_var_value ,loop_string)
% This function calls the nmf factorization but with different parms
% it loops over the set of parms provided in 
%  loop_over_var_name
%     and
%  loop_over_var_value
% it returns the scores in a nested cell array.

if isempty(loop_over_var_value)
    subsample_repeats = parms.subsample_repeats;
    scores = nan(subsample_repeats,1);
    set_terminal_title(loop_string);
    for j_sr = 1:subsample_repeats
        current_parms = parms; % to activiate parfor
        current_parms.subsample_iter = j_sr;
        samples_selected = rand_subset(j_sr,:);
        curr_X = X(samples_selected,:);
        [W, H, diff_record, time_record, eucl_dist] = load_nmf_results(curr_X, parms.num_types, current_parms.nmf_method, current_parms);

        % Match profiles to ground truth
        [W, H, best_score] = match_profiles_to_gt(W,H, GT_profiles);
        fprintf('Best mean corr is %g\n', best_score);    

        scores(j_sr) = best_score;
    end
       
else
    var_name = loop_over_var_name{1};
    loop_over_var_name = loop_over_var_name(2:end);
    var_values = loop_over_var_value{1};
    loop_over_var_value = loop_over_var_value(2:end);

    scores = {};
    for i_vars =1: length(var_values)
       current_parms = parms;
       if iscellstr(var_values)
           current_parms.(var_name) = var_values{i_vars};
           new_loop_string = sprintf('%s - %s %s',loop_string, var_name, current_parms.(var_name) );
       else
           current_parms.(var_name) = var_values(i_vars);
           new_loop_string = sprintf('%s - %s %g',loop_string, var_name, current_parms.(var_name) );
       end
       
        if strcmp(var_name, 'num_samples')
            curr_rand_subset = rand_subset{i_vars};
        else
            curr_rand_subset = rand_subset;
        end
      
        scores{i_vars} = loopOverHyperparms(X,GT_profiles,curr_rand_subset,current_parms, loop_over_var_name, loop_over_var_value ,new_loop_string); 
    end
    
end

set_terminal_title('done');
end