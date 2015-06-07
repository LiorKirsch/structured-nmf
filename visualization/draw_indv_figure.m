function draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms,y_label)
%
    FigHandle = figure('Name',parms.dataset_file); clf; hold on;
    switch getenv('USER')
      case 'lior', set(FigHandle, 'Position', [100, 100, 1600, 1000]);
    end
    set(groot,'defaultAxesColorOrder','default');       
    
    n = length(loop_over_var_value{1});    
    check_argins(loop_over_var_value, results);    
    legend_strings = loop_over_var_value{1};
    x = loop_over_var_value{2};    

    if take_from_struct(parms, 'draw_log_scale')
       x(x == 0) = 10^-5;
       x(isinf(x)) = 10^+5;
    end

    plot_h = nan(n,1);
    for i = 1:n
        curr_results = results{i};
        for j=1:3
            assert(length(curr_results) ==  length(x),'problem with build the results');
            y = cellfun(@(x) x.celltype_region_avg_scores(j), curr_results);
            subplot(1,3,j);hold on;
            plot_h(i) = plot(x, y, 's-');
        end
    end
    
    subplot(1,3,1); title('neuro');
    subplot(1,3,2); title('astro');
    subplot(1,3,3); title('oligo');
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%s %g', ...
                                               loop_over_var_name{1},x), ...
                                  legend_strings, 'UniformOutput', false);
    end
    

%     set(groot,'defaultAxesColorOrder','default');       
%     % draw the mean baselines
%     for i = 1:n
%         curr_results = results{i};
%          for j=1:3
%             subplot(1,3,j);hold on;
%             baseline = cellfun(@(x) x.baseline_celltype_region_avg_scores(j), curr_results);
%             plot_h(n +2) = plot(x, baseline,'.');
%         end
%         
%     end
%     legend_strings{end +1} = sprintf('mean-profile baseline');
%     

    set(groot,'defaultAxesColorOrder','default');       
    % draw the rand baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(1,3,j);hold on;
            set(gca,'XScale','log')
            baseline = cellfun(@(x) ...
                               x.randbase_celltype_region_avg_scores(j), ...
                               curr_results);
            baseline_sem = cellfun(@(x) ...
                                   x.randbase_celltype_region_avg_scores_sem(j), ...
                                   curr_results);
            plot_h(n +1) = errorbar(x,baseline,baseline_sem,'--');
        end
    end
    legend_strings{end +1} = sprintf('matched samples');
    legend_strings = strrep(legend_strings, '_', ' ');
    legend(plot_h, legend_strings, 'location', 'best');
    
    subplot(1,3,1);
    ylabel(y_label);
    
    for j=1:3
        subplot(1,3,j); hold on;
        
        x_label = strrep(loop_over_var_name{end}, '_', ' ') ;
        x_label = strrep(x_label, 'H lambda', '\lambda') ;
        xlabel(x_label);
        '
        switch x_label
            case '\lambda'
                xlim([10^-5.5, 10^5.5]);
                if parms.draw_log_scale 
                    skip_array = [1,2,4,6,8,9];
                    x_ticks = x(skip_array);
                    set(gca,'XTick', x_ticks);
                    tmp = get(gca,'XTickLabel');
                    tmp = strrep(tmp, '10^{-5}','indv') ;
                    tmp = strrep(tmp, '10^{5}','joined') ;
                    set(gca,'XTickLabel', tmp);
                end
            case 'num samples'

        end
%             if parms.draw_log_scale
% 
%                
%                 set(gca,'XScale','log');
%                 xticks = 10.^[-3:1:3];
%                 set(gca,'XTick', xticks);
%                 tmp = get(gca, 'XTickLabel');
%                 if any(loop_over_var_value{end} == 0 )
%                     tmp{loop_over_var_value{end} == 0} = '0';
%                 end
%                 if any(isinf(loop_over_var_value{end}))
%                     tmp{isinf(loop_over_var_value{end})} = 'inf';
%                 end
%                 set(gca,'XTickLabel', log10(xticks));
%             end
    end

    % suptitle(strrep(set_parmstr(parms), '_', ' '));
%     header = sprintf('%s-', parms.regions_short{1:end});
%     header = sprintf('%s %s:%s', header(1:end-1), parms.init_type, parms.init_subtype);
%     header = sprintf('%s NN%d', header, parms.num_restarts);
%     suptitle(header);
    parms.fig_x_axis = strrep(loop_over_var_name{end}, ' ', '_');    
    [~, file_name, dir_name] = set_filenames('figure_real', parms);
    fprintf('save figure [%s]\n', fullfile(dir_name,file_name));
    saveas(gcf,fullfile(dir_name,['indv_', file_name]));
end



function check_argins(loop_over_var_value, results);    
    recursion_levels = length(loop_over_var_value);
    msg = 'results should have 2 levels';
    assert(recursion_levels == 2 , msg);
    legend_strings = loop_over_var_value{1};
    msg = 'problem with build the results';
    assert(length(results) == length(legend_strings), msg);
end