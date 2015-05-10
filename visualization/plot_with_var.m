function plot_with_var(loop_over_var_name, loop_over_var_value, scores,score_proportions, X,GT_profiles, GT_proportions, parms)
% scores as a matrix of size [num_sizes x repeats]
% sample_sizes as a vector of size [num_sizes x 1]

   
    subplot(1,2,1);
    alg_name = draw_the_main_figure(loop_over_var_name, loop_over_var_value, scores, parms);
   
    subplot(1,2,2);
    draw_the_main_figure(loop_over_var_name, loop_over_var_value, score_proportions, parms);
    

%     switch loop_over_var_name{end}
%         case 'num_types'
%             subplot(1,2,1);hold on;
%             disp('computing random basline');
%             [baseline_scores, baseline_std, baseline_prop_score, baseline_prop_std] = get_baseline_rand_samples(X, GT_profiles,GT_proportions,loop_over_var_value{end}, parms);
%             errorbar(loop_over_var_value{end},baseline_scores,baseline_std,'k--','LineWidth',2);
%             alg_name{end +1} = 'baseline - rand sample';
% 
%             [baseline_meanProfile_scores, baseline_meanProfile_std, baseline_Mean_prop_score, baseline_Mean_prop_std] = get_baseline_profile_mean(X, GT_profiles,GT_proportions,loop_over_var_value{end},parms);
%             errorbar(loop_over_var_value{end},baseline_meanProfile_scores,baseline_meanProfile_std,'b--','LineWidth',2);
%             alg_name{end +1} = 'baseline - mean profile';
%             
%             subplot(1,2,2);hold on;
%             errorbar(loop_over_var_value{end},baseline_prop_score,baseline_prop_std,'k--','LineWidth',2);
%             errorbar(loop_over_var_value{end},baseline_Mean_prop_score,baseline_Mean_prop_std,'b--','LineWidth',2);
%         case 'num_samples'
%             disp('computing random basline');
%             num_types = loop_over_var_value{strcmp(loop_over_var_name,'num_types') };
%             
%             subplot(1,2,1);hold on;
%             [baseline_scores, baseline_std, baseline_prop_score, baseline_prop_std] = get_baseline_rand_samples(X, GT_profiles,GT_proportions,num_types, parms);
%             [baseline_scores, baseline_std, baseline_prop_score, baseline_prop_std] = do_repmat(size(loop_over_var_value{end}), baseline_scores, baseline_std, baseline_prop_score, baseline_prop_std);
%             errorbar(loop_over_var_value{end},baseline_scores,baseline_std,'k--','LineWidth',2);
%             alg_name{end +1} = 'baseline - rand sample';
% 
%             [baseline_meanProfile_scores, baseline_meanProfile_std, baseline_Mean_prop_score, baseline_Mean_prop_std] = get_baseline_profile_mean(X, GT_profiles,GT_proportions,num_types,parms);
%             [baseline_meanProfile_scores, baseline_meanProfile_std, baseline_Mean_prop_score, baseline_Mean_prop_std] = do_repmat(size(loop_over_var_value{end}), baseline_meanProfile_scores, baseline_meanProfile_std, baseline_Mean_prop_score, baseline_Mean_prop_std);
%             errorbar(loop_over_var_value{end},baseline_meanProfile_scores,baseline_meanProfile_std,'b--','LineWidth',2);
%             alg_name{end +1} = 'baseline - mean profile';
%             
%             subplot(1,2,2);hold on;
%             errorbar(loop_over_var_value{end},baseline_prop_score,baseline_prop_std,'k--','LineWidth',2);
%             errorbar(loop_over_var_value{end},baseline_Mean_prop_score,baseline_Mean_prop_std,'b--','LineWidth',2);
%     end
        
    subplot(1,2,1);        
    legend(alg_name,'Location','northwest');
    hold off;
    

    if isfield(parms,'dataset_file')
        if strncmp('okaty', parms.dataset_file,length('okaty'))
            ylim([0.8 1]);
        end
        if strncmp('barres', parms.dataset_file,length('barres'))
            switch parms.corr_type 
                case 'Spearman'
                    ylim([0.91 0.94]);
                case 'Pearson'
                    ylim([0.35 0.8]);
            end
        end
    end
    
    
    
end

function varargout = do_repmat(new_size, varargin)
    for i = 1:length(varargin)
        varargout{i} = repmat(varargin{i}, new_size );
    end
end

function alg_name = draw_the_main_figure(loop_over_var_name, loop_over_var_value, scores, parms)
    hold on;
    alg_name = loop_over_var_value{1};
    
%     assert( strcmp('nmf_method', loop_over_var_name{1} ),'nmf_method should be the first');
    
    x = loop_over_var_value{end};
    
    for i = 1:length(alg_name)
        %changing scores to - [different_parm X repeats]
        curr_scores = scores{i};
        for j=2: (length(loop_over_var_value) -1)
           curr_scores = curr_scores{1}; 
        end

        % turn the last value to the last cell
        curr_scores = cell2mat(curr_scores);
        curr_scores = curr_scores';
        
        
        y = mean(curr_scores,2);
        y_std = std(curr_scores,1,2);
        ax = errorbar(x,y,y_std,'LineWidth',2);
        
        
    end
    hold off;
    
    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    xlabel(x_label);
    ylabel(parms.corr_type);
    
    if parms.draw_log_scale
        set(gca,'XScale','log')
    end
    
    

end