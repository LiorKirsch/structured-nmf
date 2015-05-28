function draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms,y_label)
%
    FigHandle = figure('Name',parms.dataset_file); clf; hold on;
    set(FigHandle, 'Position', [100, 100, 1600, 1000]);
    set(groot,'defaultAxesColorOrder','default');       
    
    n = length(loop_over_var_value{1});
    recursion_levels = length(loop_over_var_value);
    
    assert(recursion_levels ==2 , ...
           'results should have 2 levels - the first for the legend, the second for the x axis');
    
    legend_strings = loop_over_var_value{1};
    x = loop_over_var_value{2};
    
    assert(length(results) == length(legend_strings), 'problem with build the results');

    if parms.draw_log_scale
       x(x == 0) = 10^-7;
       x(isinf(x)) = 10^+7;
    end

    plot_h = nan(n,1);
    for i = 1:n
        curr_results = results{i};
        
        for j=1:3
            assert(length(curr_results) ==  length(x),'problem with build the results');
            y = cellfun(@(x) x.celltype_region_avg_scores(j), curr_results);
            subplot(1,3,j);hold on;
            plot_h(i) = plot(x,y);
        end
    end
    
    subplot(1,3,1);
    title('neuro');
    
    subplot(1,3,2);
    title('astro');
    
    subplot(1,3,3);
    title('oligo');
    
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%s %g',loop_over_var_name{1},x), legend_strings, ...
                                  'UniformOutput', false);
    end
    

    set(groot,'defaultAxesColorOrder','default');       
    % draw the mean baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(1,3,j);hold on;
            baseline = cellfun(@(x) x.baseline_celltype_region_avg_scores(j), curr_results);
            plot_h(n +1) = plot(x,baseline,'--');
        end
        
    end
    legend_strings{end +1} = sprintf('mean-profile baseline');
    

    set(groot,'defaultAxesColorOrder','default');       
    % draw the rand baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(1,3,j);hold on;
            baseline = cellfun(@(x) x.randbase_celltype_region_avg_scores(j), curr_results);
            plot_h(n +2) = plot(x,baseline,'.');
        end
        
    end
    legend_strings{end +1} = sprintf('rand-sample baseline');

    
    legend_strings = strrep(legend_strings, '_', ' ');
    legend(plot_h, legend_strings, 'location', 'best');
    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    xlabel(x_label);
    ylabel(y_label);
    
    if parms.draw_log_scale
        
        for j=1:3
            subplot(1,3,j);hold on;
           
            set(gca,'XScale','log')
            set(gca,'XTick', x);
            tmp = get(gca,'XTickLabel');
            if any(loop_over_var_value{end} == 0 )
                tmp{loop_over_var_value{end} == 0} = '0';
            end
            if any(isinf(loop_over_var_value{end}))
                tmp{isinf(loop_over_var_value{end})} = 'inf';
            end
            set(gca,'XTickLabel', tmp);
        end
    end
    
    suptitle(strrep(set_parmstr(parms),'_',' '));
    
    parms.fig_x_axis = x_label;
    [~,file_name,dir_name] = set_filenames('figure_real', parms);
    fprintf('drawing figure - %s\n',fullfile(dir_name,file_name));

    saveas(gcf,fullfile(dir_name,['indv_', file_name]));

end