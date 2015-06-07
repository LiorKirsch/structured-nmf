function draw_figure(loop_over_var_name, loop_over_var_value, results, ...
                     parms, y_label)
%
    FigHandle = figure('Name',parms.dataset_file); clf; hold on;
    switch getenv('USER')
      case 'lior', set(FigHandle, 'Position', [100, 100, 1600, 1000]);
    end
    set(groot, 'defaultAxesColorOrder','default');       
        
    recursion_levels = length(loop_over_var_value);
    n = length(loop_over_var_value{1});

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
        assert(length(curr_results) ==  length(x),'problem with build the results');
        y = cellfun(@(x) x.run_score, curr_results);
        plot_h(i) = plot(x,y);
    end
    
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%s %g',loop_over_var_name{1},x), legend_strings, ...
                                  'UniformOutput', false);
    end

    set(groot, 'defaultAxesColorOrder','default');       
    % draw the mean baselines
    for i = 1:n
        curr_results = results{i};
        baseline = cellfun(@(x) x.baseline_score, curr_results);
        plot_h(n+1) = plot(x,baseline,'--');
    end
    legend_strings{end +1} = sprintf('mean-profile baseline');

    
    set(groot, 'defaultAxesColorOrder','default');       
    % draw the rand baselines
    for i = 1:n
        curr_results = results{i};
        baseline = cellfun(@(x) x.randbase_score, curr_results);
        baseline_sem = cellfun(@(x) x.randbase_score_sem, curr_results);
        plot_h(n+2) = errorbar(x,baseline,baseline_sem,'.');
    end
    legend_strings{end +1} = sprintf('rand-sample baseline');
        
    
    plot_h(n+3) = plot(x,repmat(0.877,size(x)),'k');
    legend_strings{end +1} = sprintf('Cortex - upper limit');
    
    legend_strings = strrep(legend_strings, '_', ' ');
    legend(plot_h,legend_strings, 'location', 'best');
    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    xlabel(x_label);
    ylabel(y_label);
    
    if parms.draw_log_scale
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
    
    title(strrep(set_parmstr(parms),'_',' '));
    
    parms.fig_x_axis = strrep(x_label, ' ', '_');
    [~, file_name, dir_name] = set_filenames('figure_real', parms);
    fprintf('save figure [%s]\n',fullfile(dir_name, file_name));
    saveas(gcf,fullfile(dir_name,file_name));
end