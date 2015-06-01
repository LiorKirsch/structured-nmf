function plot_with_var(loop_over_var_name, loop_over_var_value, scores,score_proportions, X,GT_profiles, GT_proportions, parms)
% scores as a matrix of size [num_sizes x repeats]
% sample_sizes as a vector of size [num_sizes x 1]

   
%     subplot(1,2,1);
    do_mean_corr = true;
    alg_name = draw_the_main_figure(loop_over_var_name, loop_over_var_value, scores, parms,do_mean_corr,sprintf('1 - %s',parms.corr_type));
   
%     subplot(1,2,2);
%     do_mean_corr = false;
%     draw_the_main_figure(loop_over_var_name, loop_over_var_value, score_proportions, parms,do_mean_corr,'DKL'   );
%     

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
        
%     subplot(1,2,2); 
    if ~iscell(alg_name)
        alg_name = arrayfun(@num2str, alg_name, 'UniformOutput', false);
    end
    legend(alg_name,'Location','best'); legend('boxoff');
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

function alg_name = draw_the_main_figure(loop_over_var_name, loop_over_var_value, scores, parms,do_mean_corr,y_label)
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
        
        
        if do_mean_corr
            y = nan(size(curr_scores,1),1);
            y_std = nan(size(curr_scores,1),1);
            for m =1:size(curr_scores,1);
                [y(m),y_std(m)] = mean_corr_coeff(curr_scores(m,:));
                y(m) = 1 - y(m);
            end
        else
            y = mean(curr_scores,2);
            y_std = std(curr_scores,1,2);
        end
        y_sem = y_std / sqrt(size(curr_scores,2));
        
        x(x == 0) = 10^-7;
        x(isinf(x)) = 10^+7;

        x_noise = (2*rand(size(x))-1).*x *0.05;
%         x_noise = ( i-1 )*0.5;
        ax = errorbar(x + x_noise,y,y_sem,'LineWidth',2);
        
        
    end
    hold off;

    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    xlabel(x_label);
    ylabel(y_label);
    
    if parms.draw_log_scale
        set(gca,'XScale','log')
        

        set(gca,'XTick', x);
        tmp = get(gca,'XTickLabel');
        if any(loop_over_var_value{end} == 0 )
            tmp{loop_over_var_value{end} == 0} = 'seperated models';
        end
        if any(isinf(loop_over_var_value{end}))
            tmp{isinf(loop_over_var_value{end})} = 'no region information';
        end
        set(gca,'XTickLabel', tmp);
    end
    
        

end