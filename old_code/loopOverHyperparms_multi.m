function [scores,proportions_scores] = loopOverHyperparms_multi(X, ...
                                     GT_profiles, GT_proportions,...
                                     rand_subset,parms, loop_over_var_name, ...
                                     loop_over_var_value ,loop_string)
% This function calls the nmf factorization but with different parms
% it loops over the set of parms provided in 
%  loop_over_var_name
%     and
%  loop_over_var_value
% it returns the scores in a nested cell array.

if isempty(loop_over_var_value)
    % The stopping phase of the recursion
    subsample_repeats = parms.subsample_repeats;
    scores = nan(subsample_repeats,1);
    proportions_scores = nan(subsample_repeats,1);
    set_terminal_title(loop_string);
    for j_sr = 1:subsample_repeats
        
        current_parms = parms; % to activiate parfor
        current_parms.subsample_iter = j_sr;
        
         current_parms = get_current_markers(GT_profiles,current_parms);

         
        if iscell(X)
            curr_X = cell(size(X));
            curr_GT_proportions = cell(size(X));
            for i_region = 1:length(X)
                samples_selected = rand_subset{i_region}(j_sr,:);
                curr_X{i_region} = X{i_region}(samples_selected,:);
                curr_GT_proportions{i_region} = GT_proportions{i_region}(:,samples_selected);
            end
            
        else
            samples_selected = rand_subset(j_sr,:);
            curr_X = X(samples_selected,:);
            curr_GT_proportions = GT_proportions(:,samples_selected);
        end
        
        [W, H, diff_record, time_record, eucl_dist] = ...
            load_nmf_results(curr_X, parms.num_types, ...
                             current_parms.nmf_method, current_parms);
                         
        if iscell(X)
            best_scores = nan(length(X),1);
            best_proportions_scores = nan(length(X),1);
            for i_region = 1:length(X)
                % Match profiles to ground truth
                %TOTDO ADD LOOP
               
                [W{i_region}, H{i_region}, best_scores(i_region), best_proportions_scores(i_region)] = match_profiles_to_gt(...
                    W{i_region},H{i_region}, ...
                    GT_profiles{i_region}, curr_GT_proportions{i_region}, ...
                    parms.corr_type);
                fprintf('Best mean corr (%s) is %g (proprtions %g)\n', ...
                    parms.regions{i_region},best_scores(i_region),best_proportions_scores(i_region));    

            end
            best_score = mean_corr_coeff(best_scores);
            best_proportions_score = mean(best_proportions_scores);
            
        else
            % Match profiles to ground truth
            [W, H, best_score, best_proportions_score] = match_profiles_to_gt(W,H, ...
                GT_profiles, curr_GT_proportions, parms.corr_type);
            fprintf('Best mean corr is %g (proprtions %g)\n', ...
                best_score,best_proportions_score);    
        
        end
        
        scores(j_sr) = best_score;
        proportions_scores(j_sr) = best_proportions_score;
        
    end

%     figure;
%     ax1 = subplot(1,2,1);
%     draw_profiles(GT_profiles,GT_proportions,parms);
%     ax2 = subplot(1,2,2);
%     draw_profiles_with_GT(H,cellfun(@transpose ,W,'UniformOutput',false),GT_profiles,parms);
%     linkaxes([ax1,ax2],'xy');

%    figure; draw_profiles(GT_profiles,cellfun(@transpose ,W,'UniformOutput',false),parms);
%    figure; draw_profiles(H,cellfun(@transpose ,W,'UniformOutput',false),parms);
%     title(set_parmstr(parms));
else
    % Remove one layer from the recursion
    var_name = loop_over_var_name{1};
    loop_over_var_name = loop_over_var_name(2:end);
    var_values = loop_over_var_value{1};
    loop_over_var_value = loop_over_var_value(2:end);

    scores = {};
    proportions_scores = {};
    for i_vars = 1: length(var_values)
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
      
        [scores{i_vars}, proportions_scores{i_vars}] = loopOverHyperparms_multi(X,GT_profiles, ...
                                            GT_proportions, ...
                                            curr_rand_subset,current_parms, ...
                                            loop_over_var_name, ...
                                            loop_over_var_value, ...
                                            new_loop_string); 
    end
    
end

set_terminal_title('done');
end

function parms = get_current_markers(GT_profiles,parms)

     if isfield(parms,'num_markers')
         
         if iscell(GT_profiles)
             
            H_markers = cellfun(@(x) get_markers_from_GT(x', parms.num_markers, 1000,parms),...
                GT_profiles,'UniformOutput',false);
            parms.H_markers = H_markers;  
            
        else
            H_markers = get_markers_from_GT(GT_profiles', parms.num_markers, 1000,parms);
            parms.H_markers = H_markers;
         end
        
    end
    
end